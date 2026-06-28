#!/usr/bin/env python3
"""Convert common kinetics expressions to VULCAN A, B, C parameters."""

from __future__ import annotations

import argparse
import math

AVOGADRO = 6.02214076e23
R_BY_UNIT = {
    "j/mol": 8.31446261815324,
    "kj/mol": 8.31446261815324e-3,
    "cal/mol": 1.98720425864083,
    "kcal/mol": 1.98720425864083e-3,
    "k": 1.0,
}


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--temperatures",
        nargs="*",
        type=float,
        default=[],
        help="Temperatures at which to evaluate the converted rate",
    )


def molecule_a(a_value: float, order: int, quantity_basis: str) -> float:
    if order < 1:
        raise ValueError("Reaction order must be at least 1")
    if quantity_basis == "molecule":
        return a_value
    return a_value / AVOGADRO ** (order - 1)


def activation_temperature(energy: float, unit: str) -> float:
    return energy / R_BY_UNIT[unit]


def evaluate(a_value: float, b_value: float, c_value: float, temperature: float) -> float:
    if temperature <= 0:
        raise ValueError("Temperature must be positive")
    return a_value * temperature**b_value * math.exp(-c_value / temperature)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="kind", required=True)

    nist = sub.add_parser("nist", help="Convert Aref*(T/Tref)^n*exp(-Ea/RT)")
    nist.add_argument("--a-ref", type=float, required=True)
    nist.add_argument("--t-ref", type=float, default=298.0)
    nist.add_argument("--n", type=float, required=True)
    nist.add_argument("--ea", type=float, required=True)
    nist.add_argument("--ea-unit", choices=R_BY_UNIT, default="j/mol")
    nist.add_argument("--order", type=int, default=2)
    nist.add_argument("--quantity-basis", choices=("molecule", "mol"), default="molecule")
    add_common(nist)

    kooij = sub.add_parser("kooij", help="Convert alpha*(T/Tref)^beta*exp(-gamma/T)")
    kooij.add_argument("--alpha", type=float, required=True)
    kooij.add_argument("--t-ref", type=float, default=300.0)
    kooij.add_argument("--beta", type=float, required=True)
    kooij.add_argument("--gamma", type=float, required=True)
    kooij.add_argument("--order", type=int, default=2)
    kooij.add_argument("--quantity-basis", choices=("molecule", "mol"), default="molecule")
    add_common(kooij)

    chemkin = sub.add_parser("chemkin", help="Convert A*T^b*exp(-Ea/RT)")
    chemkin.add_argument("--a", type=float, required=True)
    chemkin.add_argument("--b", type=float, required=True)
    chemkin.add_argument("--ea", type=float, required=True)
    chemkin.add_argument("--ea-unit", choices=R_BY_UNIT, default="cal/mol")
    chemkin.add_argument("--order", type=int, required=True)
    chemkin.add_argument("--quantity-basis", choices=("molecule", "mol"), default="mol")
    add_common(chemkin)

    return parser


def main() -> None:
    args = build_parser().parse_args()

    if args.kind == "nist":
        a_value = molecule_a(args.a_ref, args.order, args.quantity_basis) / args.t_ref**args.n
        b_value = args.n
        c_value = activation_temperature(args.ea, args.ea_unit)
    elif args.kind == "kooij":
        a_value = molecule_a(args.alpha, args.order, args.quantity_basis) / args.t_ref**args.beta
        b_value = args.beta
        c_value = args.gamma
    else:
        a_value = molecule_a(args.a, args.order, args.quantity_basis)
        b_value = args.b
        c_value = activation_temperature(args.ea, args.ea_unit)

    print(f"A = {a_value:.10E}")
    print(f"B = {b_value:.10g}")
    print(f"C = {c_value:.10g} K")
    print(f"VULCAN columns: {a_value:.8E} {b_value:.8g} {c_value:.8g}")
    for temperature in args.temperatures:
        value = evaluate(a_value, b_value, c_value, temperature)
        print(f"k({temperature:g} K) = {value:.10E}")


if __name__ == "__main__":
    main()

