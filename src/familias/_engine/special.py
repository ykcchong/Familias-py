"""Helpers from C++ ``special.cpp`` we still need in Python."""
from __future__ import annotations

# Maximum entries in a cutset memoisation table (port of MAX_TABLE_SIZE).
MAX_TABLE_SIZE: int = 100000


def get_name_prefix(names) -> str:
    """Return a string prefix such that ``prefix + str(int)`` does not
    collide with any name in ``names``.

    Faithful behaviour of ``special.cpp::getNamePrefix``: starts with "E"
    and prepends additional "E" characters whenever the prefix is itself
    a prefix of one of the names.
    """
    prefix = "E"
    while True:
        if not any(n.startswith(prefix) for n in names):
            return prefix
        prefix = "E" + prefix


def mypow(x: float, n: int) -> float:
    """Integer-exponent power, identical semantics to ``special::mypow``."""
    if n == 0:
        return 1.0
    if n < 0:
        return 1.0 / mypow(x, -n)
    out = 1.0
    for _ in range(n):
        out *= x
    return out
