from pathlib import Path
import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation
import anndata as ad
import gc

# ------------------------------
# 你现有的工具函数（略作健壮性优化）
# ------------------------------
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def find_matrix_dir(sample_dir: Path) -> Path:
    candidates = [
        sample_dir / "outs/Gene/filtered",                 # STARsolo
        sample_dir / "outs/filtered_feature_bc_matrix",    # Cell Ranger (filtered)
    ]
    for p in candidates:
        if p.is_dir():
            return p
    raise FileNotFoundError(f"No 10x matrix dir under {sample_dir}")

def read_one_sample(sample_id: str, BASE: Path) -> ad.AnnData:
    sdir = BASE / sample_id
    mdir = find_matrix_dir(sdir)

    # 注意：10x矩阵常是barcode在obs，基因在var
    adata = sc.read_10x_mtx(mdir, var_names='gene_symbols')
    adata.var_names_make_unique()

    # obs 命名与标签
    adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"]  = sample_id
    adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs["barcode"]]

    # QC 变量
    g = adata.var_names.str.upper()
    adata.var["mt"]   = g.str.startswith("MT-")
    adata.var["ribo"] = g.str.startswith(("RPS", "RPL"))
    adata.var["hb"]   = g.str.match(r"^HB(?!P)")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, inplace=True
    )
    return adata

# ------------------------------
# 导出与分批合并主逻辑
# ------------------------------
def export_each_sample_h5ad(sample_ids, BASE: Path, OUTDIR: Path, overwrite: bool=False):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for sid in sample_ids:
        out_path = OUTDIR / f"{sid}.h5ad"
        if out_path.exists() and not overwrite:
            print(f"[skip] {sid} 已存在：{out_path}")
            continue
        print(f"[read/write] 读取 {sid} 并写出 {out_path}")
        adata = read_one_sample(sid, BASE)
        # 推荐开启压缩，节省空间；也可不压缩以加快速度
        adata.write(out_path, compression="gzip")
        del adata
        gc.collect()

def chunk_list(lst, n_chunks: int):
    # 将列表均匀分成n_chunks份
    return np.array_split(list(lst), n_chunks)

def concat_h5ad_files(h5ad_paths, label: str, keys, out_path: Path):
    """将一批h5ad文件读入合并并写出"""
    adatas = []
    for p in h5ad_paths:
        adata = ad.read_h5ad(p)
        adatas.append(adata)
    print(f"[concat] 合并 {len(adatas)} 个对象 → {out_path.name}")
    adata_batch = ad.concat(
        adatas, join="outer", label=label, keys=keys, index_unique=None
    )
    # 写出并释放内存
    adata_batch.write(out_path, compression="gzip")
    for a in adatas:
        a.file.close() if a.isbacked else None
        del a
    del adata_batch
    gc.collect()

def batched_merge(sample_ids, PER_SAMPLE_DIR: Path, BATCH_DIR: Path, FINAL_PATH: Path, n_batches: int=4, overwrite: bool=False):
    """
    1) 将 {PER_SAMPLE_DIR}/{sample}.h5ad 分成 n_batches 组，分别合并为 batch_i.h5ad
    2) 最终再将 n_batches 个 batch 合并为 FINAL_PATH
    """
    BATCH_DIR.mkdir(parents=True, exist_ok=True)
    chunks = chunk_list(sample_ids, n_batches)

    batch_paths = []
    for i, chunk in enumerate(chunks, start=1):
        chunk = list(map(str, chunk))  # numpy array -> list[str]
        out_path = BATCH_DIR / f"batch_{i}.h5ad"
        batch_paths.append(out_path)
        if out_path.exists() and not overwrite:
            print(f"[skip] 批次 {i} 已存在：{out_path}")
            continue

        # 该批次的单样本文件路径
        h5ads = [PER_SAMPLE_DIR / f"{sid}.h5ad" for sid in chunk]
        missing = [str(p) for p in h5ads if not p.exists()]
        if missing:
            raise FileNotFoundError(f"以下单样本文件缺失，先运行导出：\n" + "\n".join(missing))

        # 合并该批次
        concat_h5ad_files(
            h5ad_paths=h5ads,
            label="sample",
            keys=chunk,
            out_path=out_path
        )

    # 最终合并四个批次
    if FINAL_PATH.exists() and not overwrite:
        print(f"[skip] 最终合并文件已存在：{FINAL_PATH}")
        return

    print(f"[final concat] 合并 {len(batch_paths)} 个批次 → {FINAL_PATH.name}")
    concat_h5ad_files(
        h5ad_paths=batch_paths,
        label="batch",
        keys=[p.stem for p in batch_paths],  # e.g., batch_1, batch_2, ...
        out_path=FINAL_PATH
    )

# ------------------------------
# 参数与执行
# ------------------------------
BASE = Path("/home/sduan/scrup/alignments")
PER_SAMPLE_DIR = Path("/home/sduan/scrup/object/scanpy")     # 1) 每个样本导出目录
BATCH_DIR      = Path("/home/sduan/scrup/object/scanpy_batches")  # 2) 批次合并目录
FINAL_PATH     = Path("/home/sduan/scrup/object/scanpy_all.h5ad") # 3) 最终合并结果

# 枚举样本（SRX开头目录）
sample_ids = sorted([p.name for p in BASE.iterdir() if p.is_dir() and p.name.startswith("SRX")])

# 先逐样本读取并写出单独的 h5ad（可多次运行，已存在则跳过）
export_each_sample_h5ad(sample_ids, BASE, PER_SAMPLE_DIR, overwrite=False)

# 将900+样本分为4个批次合并，再最终合并四个批次
batched_merge(sample_ids, PER_SAMPLE_DIR, BATCH_DIR, FINAL_PATH, n_batches=4, overwrite=False)

print("[done] 单样本导出、4批次合并以及最终合并完成。")
