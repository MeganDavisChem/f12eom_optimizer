"""Console script for f12eom_optimizer."""
import sys
import click


@click.command()
def main(args=None):
    import f12eom_optimizer.f12eom_optimizer as f12
    f12.main()
    """Console script for f12eom_optimizer."""
    click.echo("Replace this message by putting your code into "
               "f12eom_optimizer.cli.main")
    click.echo("See click documentation at https://click.palletsprojects.com/")
    return 0


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
