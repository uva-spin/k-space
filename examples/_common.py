"""Shared validation and reproducibility helpers for runnable scientific examples."""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
# Run directly from a source checkout; no editable install is required.
sys.path.insert(0, str(ROOT / 'src'))
from resummed_operator import CanonicalRG, Coefficients  # noqa: E402


def finite_float(value: str) -> float:
    import math
    number = float(value)
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError('A finite number is required')
    return number


def add_rg_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument('--coefficients', type=Path, default=ROOT/'data/coefficients_demo.json')
    parser.add_argument('--mu', type=finite_float, default=20.0, help='Final virtuality scale in GeV')
    parser.add_argument('--zeta', type=finite_float, default=400.0, help='Final rapidity scale in GeV^2')
    parser.add_argument('--mu0', type=finite_float, default=2.0, help='Infrared completion scale in GeV')
    parser.add_argument('--a0', type=finite_float, default=0.02, help='alpha_s(mu0)/(4*pi), not alpha_s itself')


def build_rg(args: argparse.Namespace) -> CanonicalRG:
    if args.mu0 <= 0 or args.a0 <= 0 or args.mu < args.mu0 or args.zeta <= 0:
        raise ValueError('Require mu0>0, a0>0, mu>=mu0, and zeta>0')
    coefficients = Coefficients.load(args.coefficients)
    if coefficients.physical_n3llprime_payload:
        raise ValueError('These examples use the documented demonstration input; metadata alone cannot authorize a physical N3LL-prime result')
    return CanonicalRG(coefficients, mu0=args.mu0, a0=args.a0)


def output_directory(path: Path) -> Path:
    """Write only below this checkout's outputs tree; never delete existing files."""
    original = path.expanduser().absolute()
    if any(parent.is_symlink() for parent in [original, *original.parents]):
        raise ValueError('Symlinks are not permitted in example output paths')
    resolved = original.resolve()
    allowed = (ROOT/'outputs').resolve()
    if resolved == allowed or not resolved.is_relative_to(allowed):
        raise ValueError('Choose a subdirectory of this checkout\'s outputs/ directory')
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def provenance(args: argparse.Namespace) -> dict:
    return {
        'coefficients_filename': args.coefficients.name,
        'coefficients_sha256': hashlib.sha256(args.coefficients.read_bytes()).hexdigest(),
        'operator_sha256': hashlib.sha256((ROOT/'src/resummed_operator.py').read_bytes()).hexdigest(),
        'python': platform.python_version(),
        'dependencies': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy')},
        'scales': {'mu_GeV': args.mu, 'zeta_GeV2': args.zeta, 'mu0_GeV': args.mu0, 'a0': args.a0},
        'scope': {
            'coefficient_input': 'declared demonstration input; see docs/COEFFICIENTS.md',
            'physical_N3LLprime_prediction_validated': False,
            'continuum_grid_error_certified': False,
            'data_fit_performed': False,
        },
    }


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8')
