"""Build the database for dram2"""
from pathlib import Path
import click
import yaml
import pkgutil
import importlib


from dram2.cli.context import (
    DramContext,
    log_error_wraper,
    get_time_stamp_id,
    __version__,
    USER_CONFIG,
)
from dram2.db_kits.utils import (
    DBKit,
)
from dram2.utils import DramUsageError

from dram2 import db_kits as db_kits


for i in pkgutil.iter_modules(db_kits.__path__, db_kits.__name__ + "."):
    importlib.import_module(i.name)

DB_KITS: list = [
    i for i in DBKit.__subclasses__() if i.selectable and len(i.location_keys) > 0
]


@click.command(
    "build_db_list",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
def db_list_cmd():
    """
    List the input files you can provied
    ___

    Just needed a whay to show the input file options for the build_db comand. Use the
    -i argument with the name and the path.

    use::

        conda activate ./dram2_env
        dram2 build_db_list
    """
    for i in DB_KITS:
        print(f"db {i}:")
        for j in i.location_keys:
            print(f"   {j}")


@click.command(
    "build_db",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.option(
    "--make_db",
    multiple=True,
    default=[i.name for i in DB_KITS],
    type=click.Choice([i.name for i in DB_KITS]),
    help="""
    Specify exactly which DBs to use. This argument can be used multiple times, so
    for example if you want to annotate with FeGenie and Camper you would have a
    command like `dram2 - o output/dir annotate - -use_db fegenie - -use_db
    camper`, the options available are in this help.
    """,
)
@click.option(
    "-y",
    "--output_yaml",
    type=click.Path(
        path_type=Path,
    ),
    default=USER_CONFIG,
    help="""
    Location where you would like dram to put the yaml that tells dram where to
    find all the data. There are pre-set locations were DRAM looks for the yaml.
    With out this argument DRAM2 will put the config in
    `USER_HOME/.config/dram_config.yaml`, only if it is empty. If after testing the
    user config you decide that you like it`/etc/dram2_config.yaml`.

    """,
)
@click.option(
    "-u",
    "--update",
    is_flag=True,
    show_default=True,
    default=False,
    help="""
    Update any database that is requested, this will erase what is already in the
    yaml for that database. If a yaml is not specified
    then the user database is read and updated.
    If you want to update the global
    you need to point to it with the yaml argument and run with the needed
    perminsion. The global config is /etc/dram2_config.yaml.
    """,
)
@click.option(
    "-i",
    "--input_file",
    nargs=2,
    type=click.Tuple(
        [
            click.Choice(
                [j for i in DB_KITS for j in i.location_keys],
                case_sensitive=False
                # [j for i in DB_KITS for j in i.location_keys], case_sensitive=False
            ),
            click.Path(exists=True, path_type=Path),
        ]
    ),
    multiple=True,
    help="""
    Sometimes you want to provide database files your self. Maybe you have a custom
    version of the file, maybe you can't download a file or just don't want too. Use
    this argumemet to point to those files with the format "-i <name> <path>
    """,
)
@click.pass_context
def db_builder_cmd(
    ctx: click.Context,
    make_db: list[str],
    update: bool,
    input_file: tuple[str, Path],
    output_yaml: Path = USER_CONFIG,
):
    """
    Build Your Own Custom DRAM Database
    ___

    A tool to build your own custom dram database. Download the raw data and format
    it for use in DRAM2. This is only for advanced users but it is necessary for our
    kegg users.

    All your data will apear in the DRAM-dir wich you specify with the root command. If
    you want to move that folder you will need to edit the yaml to point to the new
    location, so try to pick it carefully.

    If you want to update the global
    you need to point to it with the yaml argument and run with the needed
    perminsion. The global config is /etc/dram2_config.yaml.

    Example Use::

        conda activate ./dram2_env
        pip install -e ./
        dram2 build_db -h
        dram2 -d temp4 build_db  -y test.yaml
        dram2 -d temp4 build_db  -y test.yaml --update


    TODO:

    - Fully Implement the Download options in the DB_Kits
    - Fully Implement the process options in the DB_KITS

    """
    context: DramContext = ctx.obj
    log_error_wraper(db_builder_pipe, context, "database builder")(
        context=context,
        make_db=make_db,
        update=update,
        input_file=input_file,
        output_yaml=output_yaml,
    )


def db_builder_pipe(
    context: DramContext,
    make_db: list[str],
    update: bool,
    input_file: list[tuple[str, Path]],
    output_yaml: Path,
):
    """
    Build Your Own Custom DRAM Database (Not Ready)
    ___

    A tool to build your own custom dram database. Download the raw data and format
    it for use in dram. This is only for advanced users but it is necessary for our
    kegg users.

    This command should only call process that are in the db_kits. That is where
    downloading should be handled.

    TODO:

    - Fully Implement the Download options in the DB_Kits
    - Fully Implement the process options in the DB_KITS

    """
    dram_dir: Path = context.get_dram_dir()
    logger = context.get_logger()
    databases: list[DBKit] = [i({}, logger) for i in DB_KITS if i.name in set(make_db)]
    if len(databases) < 1:
        raise DramUsageError("No data bases to build")
    input_dict: dict[str, Path] = {i: j for i, j in input_file}
    for db in databases:
        db_builder(
            db,
            input_dict,
            dram_dir,
            output_yaml,
            update,
        )


def db_builder(
    db_kit: DBKit,
    input_dict: dict[str, Path],
    out_dir: Path,
    out_yaml: Path,
    update: bool = False,
):
    if out_yaml.exists():
        with open(out_yaml, "r") as conf:
            config = yaml.safe_load(conf)
    else:
        config = {}
    if isinstance(conf_folder_str := config.get("dram_data_folder"), str):
        if Path(conf_folder_str) != out_dir.absolute():
            raise DramUsageError(
                f"The DRAM-config file has an output directory: {conf_folder_str} "
                f"but the command has the output dir {out_dir.absolute()}. This is "
                f"not suported by this script. Please put all the data in one folder "
                f"or don't use a pre-existing config."
            )
    name = db_kit.name
    db_kit.set_args(out_dir / name)
    download_paths: dict[str, Path] = {}
    download_paths.update(db_kit.download(input_dict))
    config_out = db_kit.setup(download_paths, out_dir)
    if "db_kits" not in config:
        config.update({"db_kits": {}})
    if name in config["db_kits"] and not update:
       raise DramUsageError(
            f"There is already data for {name}, and you have not "
            f"specified this as an update with the update flag so DRAM2 "
            f"will not erase that data from the config file."
        )

    config.update(
        {
            "dram_data_folder": out_dir.absolute().as_posix(),
            "dram_version": __version__,
            "log_file_path": "dram2.log",
        }
    )
    config["db_kits"].update({name: config_out})
    with open(out_yaml, "w") as dram_conf:
        yaml.safe_dump(config, dram_conf)
