"Test comandline arguments using mocker"
#! /bin/python3
import pytest

def test_arguments_from_cli(mocker):
    """Test whether arguments from the command line are set up correctly in a HermesApp object."""
    mocker.patch(
        "sys.argv",
        [
            "rhasspy-hermes-app-test",
            "--host",
            "rhasspy.home",
            "--port",
            "8883",
            "--tls",
            "--username",
            "rhasspy-hermes-app",
            "--password",
            "test",
        ],
    )
    app = HermesApp("Test arguments in init", mqtt_client=mocker.MagicMock())

    assert app.args.host == "rhasspy.home"
    assert app.args.port == 8883
    assert app.args.tls is True
    assert app.args.username == "rhasspy-hermes-app"
    assert app.args.password == "test"

