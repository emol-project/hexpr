"""hexpr status code -> Python exception mapping."""


class HexprError(RuntimeError):
    """Raised when a hexpr_* function returns non-zero."""


class InvalidArg(HexprError):
    pass


class OutOfMemory(HexprError):
    pass


class IOError(HexprError):
    pass


class RangeError(HexprError):
    pass


class NotBuilt(HexprError):
    pass


class InternalError(HexprError):
    pass


_STATUS_MAP = {
    1:  InvalidArg,
    2:  OutOfMemory,
    3:  IOError,
    4:  RangeError,
    5:  NotBuilt,
    99: InternalError,
}


def check_status(status, lib=None):
    if status == 0:
        return
    cls = _STATUS_MAP.get(status, HexprError)
    msg = f"hexpr returned status {status}"
    if lib is not None:
        raw = lib.hexpr_last_error()
        if raw:
            msg += f": {raw.decode('utf-8', errors='replace')}"
    raise cls(msg)
