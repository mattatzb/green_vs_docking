"""One entry point, with optional commands for individual pipeline stages."""

from __future__ import annotations

import argparse
import importlib
import sys

COMMANDS = {"run": "run", "score": "score", "plot": "plot", "table": "table"}


def main(argv: list[str] | None = None) -> None:
    arguments = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        description="Manuscript analysis pipeline. With no command, run the complete workflow.",
        epilog="Use 'COMMAND --help' for stage options. Generated results default to outputs/9KQ3.",
    )
    parser.add_argument("command", nargs="?", choices=COMMANDS, default="run",
                        help="run: full pipeline; score: one comparison; plot: plot/resume; table: final table")
    if arguments and arguments[0] in {"-h", "--help"}:
        parser.parse_args(arguments)
        return
    if arguments and arguments[0] in COMMANDS:
        command = arguments.pop(0)
    elif not arguments or arguments[0].startswith("-"):
        command = "run"
    else:
        parser.error(f"unknown command: {arguments[0]}")
    module = importlib.import_module(f".{COMMANDS[command]}", package=__package__)
    module.main(arguments)
