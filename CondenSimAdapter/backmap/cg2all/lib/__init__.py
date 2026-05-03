__all__ = ["convert_cg2all"]


def __getattr__(name: str):
    if name == "convert_cg2all":
        from .snippets import convert_cg2all
        return convert_cg2all
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")