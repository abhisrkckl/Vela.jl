from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

from .spnta import SPNTA


def parse_args(argv):
    parser = ArgumentParser(
        prog="pyvela-jlso",
        description="Read a par file, tim file, and prior JSON file, and write a JLSO file.",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "par_file",
        help="The pulsar ephemeris file. Should be readable using PINT. The "
        "uncertainties listed in the file will be used for 'cheat' priors where applicable.",
    )
    parser.add_argument(
        "tim_file", help="The pulsar TOA file. Should be readable using PINT."
    )
    parser.add_argument(
        "-P",
        "--prior_file",
        help="A JSON file containing the prior distributions for each free parameter.",
    )
    parser.add_argument(
        "-C",
        "--cheat_prior_scale",
        default=50,
        type=float,
        help="The scale factor by which the frequentist uncertainties are multiplied to "
        "get the 'cheat' prior distributions.",
    )
    parser.add_argument(
        "--no_marginalize_gp_noise",
        action="store_true",
        help="Don't analytically marginalize the correlated Gaussian noise amplitudes.",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="The output file name. Will replace an existing file.",
        required=True,
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    spnta = SPNTA(
        args.par_file,
        args.tim_file,
        cheat_prior_scale=args.cheat_prior_scale,
        custom_priors=(args.prior_file if args.prior_file is not None else {}),
        marginalize_gp_noise=(not args.no_marginalize_gp_noise),
    )
    spnta.save_jlso(args.outfile)
