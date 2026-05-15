"""Internal module: GROMACS topology parsing helpers for softcore optimization."""

import os
import re
import shutil

novarcharre = re.compile(r"\W")


def _find_all_instances_in_string(string, substr):
    """Find indices of all instances of substr in string"""
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx + 1)
    return indices


def _replace_defines(line, defines):
    """Replaces defined tokens in a given line"""
    if not defines:
        return line
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices:
            continue
        # Check to see if it's inside of quotes
        inside = ""
        idx = 0
        n_to_skip = 0
        new_line = []
        for i, char in enumerate(line):
            if n_to_skip:
                n_to_skip -= 1
                continue
            if char in ("'\""):
                if not inside:
                    inside = char
                else:
                    if inside == char:
                        inside = ""
            if idx < len(indices) and i == indices[idx]:
                if inside:
                    new_line.append(char)
                    idx += 1
                    continue
                if i == 0 or novarcharre.match(line[i - 1]):
                    endidx = indices[idx] + len(define)
                    if endidx >= len(line) or novarcharre.match(line[endidx]):
                        new_line.extend(list(value))
                        n_to_skip = len(define) - 1
                        idx += 1
                        continue
                idx += 1
            new_line.append(char)
        line = "".join(new_line)

    return line


def _defaultGromacsIncludeDir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top"""
    if "GMXDATA" in os.environ:
        return os.path.join(os.environ["GMXDATA"], "top")
    if "GMXBIN" in os.environ:
        return os.path.abspath(os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top"))

    pdb2gmx_path = shutil.which("pdb2gmx")
    if pdb2gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(pdb2gmx_path), "..", "share", "gromacs", "top")
        )

    gmx_path = shutil.which("gmx")
    if gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(gmx_path), "..", "share", "gromacs", "top")
        )

    return "/usr/local/gromacs/share/gromacs/top"
