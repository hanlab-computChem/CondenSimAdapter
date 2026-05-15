"""Unified subprocess runner for external command execution.

All external command invocations in the minimize package should go through
this module so that capture_output, error messages, and "command not found"
handling are consistent.
"""

import subprocess
from pathlib import Path
from typing import List, Optional


class SubprocessError(RuntimeError):
    """Raised when an external command fails."""

    def __init__(self, cmd: List[str], returncode: int, stderr: str, cwd: Optional[Path] = None):
        self.cmd = cmd
        self.returncode = returncode
        self.stderr = stderr
        self.cwd = cwd
        msg = f"Command '{cmd[0]}' failed (code {returncode})"
        if cwd:
            msg += f" in {cwd}"
        msg += f": {stderr[-500:]}"
        super().__init__(msg)


def _run_command(
    cmd: List[str],
    *,
    cwd: Optional[Path] = None,
    input_str: Optional[str] = None,
    timeout: Optional[float] = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run an external command with unified error handling.

    Args:
        cmd: Command and arguments as a list of strings.
        cwd: Working directory for the subprocess.
        input_str: Text to pipe to stdin.
        timeout: Timeout in seconds.
        check: If True, raise SubprocessError on non-zero exit.

    Returns:
        The completed process.

    Raises:
        FileNotFoundError: If the command binary is not found.
        SubprocessError: If check=True and the command exits non-zero.
    """
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            input=input_str,
            cwd=str(cwd) if cwd else None,
            timeout=timeout,
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Command '{cmd[0]}' not found. Ensure it is installed and on PATH."
        )

    if check and result.returncode != 0:
        raise SubprocessError(cmd, result.returncode, result.stderr, cwd)

    return result
