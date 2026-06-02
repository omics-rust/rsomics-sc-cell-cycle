use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const EPSILON: f64 = 1e-5;

struct Row {
    s: f64,
    g2m: f64,
    phase: String,
}

fn scanpy_python() -> Option<PathBuf> {
    let home = std::env::var("HOME").ok()?;
    let shared = PathBuf::from(&home).join("oracle-venvs/scanpy/bin/python");
    if shared.exists() {
        return Some(shared);
    }
    for cand in ["python3", "python"] {
        let ok = Command::new(cand)
            .args(["-c", "import scanpy"])
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if ok {
            return Some(PathBuf::from(cand));
        }
    }
    None
}

fn parse_tsv(text: &str) -> HashMap<String, Row> {
    let mut m = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let mut it = line.split('\t');
        let id = it.next().unwrap().to_string();
        let s = it.next().unwrap().parse::<f64>().unwrap();
        let g2m = it.next().unwrap().parse::<f64>().unwrap();
        let phase = it.next().unwrap().to_string();
        m.insert(id, Row { s, g2m, phase });
    }
    m
}

fn diff(want: &HashMap<String, Row>, got: &HashMap<String, Row>, label: &str) -> f64 {
    assert_eq!(
        want.len(),
        got.len(),
        "{label}: {} vs {} cells",
        want.len(),
        got.len()
    );
    let mut max_dev = 0.0_f64;
    for (id, w) in want {
        let g = got
            .get(id)
            .unwrap_or_else(|| panic!("{label}: missing cell {id}"));
        let ds = (w.s - g.s).abs();
        let dg = (w.g2m - g.g2m).abs();
        max_dev = max_dev.max(ds).max(dg);
        assert!(
            ds < EPSILON,
            "{label}: {id} S_score differs: {} vs {}",
            w.s,
            g.s
        );
        assert!(
            dg < EPSILON,
            "{label}: {id} G2M_score differs: {} vs {}",
            w.g2m,
            g.g2m
        );
        assert_eq!(w.phase, g.phase, "{label}: {id} phase differs");
    }
    max_dev
}

fn run_ours(mtx_dir: &Path, s: &Path, g2m: &Path, extra: &[&str]) -> HashMap<String, Row> {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rsomics-sc-cell-cycle"));
    cmd.arg(mtx_dir)
        .arg("--s-genes")
        .arg(s)
        .arg("--g2m-genes")
        .arg(g2m)
        .arg("-o")
        .arg("-")
        .arg("-q");
    cmd.args(extra);
    let out = cmd.output().expect("run rsomics-sc-cell-cycle");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    parse_tsv(&String::from_utf8(out.stdout).unwrap())
}

#[test]
fn matches_committed_golden() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let mtx_dir = Path::new(manifest).join("tests/golden/tenx");
    let s = Path::new(manifest).join("tests/golden/s_genes.txt");
    let g2m = Path::new(manifest).join("tests/golden/g2m_genes.txt");
    let golden = Path::new(manifest).join("tests/golden/cell_cycle.tsv");
    assert!(golden.exists(), "missing golden {golden:?}");

    let ours = run_ours(&mtx_dir, &s, &g2m, &[]);
    let want = parse_tsv(&std::fs::read_to_string(&golden).unwrap());
    let dev = diff(&want, &ours, "golden");
    eprintln!("golden compat OK: {} cells, max dev {dev:e}", ours.len());
}

#[test]
fn matches_scanpy_value_level() {
    let Some(py) = scanpy_python() else {
        eprintln!("SKIP: scanpy venv not found (~/oracle-venvs/scanpy/bin/python); compat skipped");
        return;
    };

    let manifest = env!("CARGO_MANIFEST_DIR");
    let mtx_dir = Path::new(manifest).join("tests/golden/tenx");
    let s = Path::new(manifest).join("tests/golden/s_genes.txt");
    let g2m = Path::new(manifest).join("tests/golden/g2m_genes.txt");
    let oracle_py = Path::new(manifest).join("tests/scanpy_cell_cycle_oracle.py");
    assert!(mtx_dir.exists(), "missing golden 10x dir {mtx_dir:?}");

    let scratch = std::env::var("RSOMICS_SCRATCH")
        .map(PathBuf::from)
        .unwrap_or_else(|_| std::env::temp_dir());

    let oracle_out = scratch.join("sc_cell_cycle_oracle_default.tsv");
    let status = Command::new(&py)
        .arg(&oracle_py)
        .arg(&mtx_dir)
        .arg(&s)
        .arg(&g2m)
        .arg(&oracle_out)
        .status()
        .expect("run scanpy oracle");
    assert!(status.success(), "oracle failed");
    let oracle = parse_tsv(&std::fs::read_to_string(&oracle_out).unwrap());
    let ours = run_ours(&mtx_dir, &s, &g2m, &[]);
    let dev = diff(&oracle, &ours, "default");
    eprintln!("live compat OK: {} cells, max dev {dev:e}", ours.len());
}
