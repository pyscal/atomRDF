from typing import List, Optional, Union, Generic, TypeVar, get_args, get_origin, Any
from pydantic import BaseModel as PydanticBaseModel
from pydantic import Field

T = TypeVar("T")


def remove_empty_dicts(data):
    if isinstance(data, dict):
        cleaned = {
            k: remove_empty_dicts(v)
            for k, v in data.items()
            if not (isinstance(v, dict) and remove_empty_dicts(v) == {})
        }
        return cleaned
    elif isinstance(data, list):
        return [remove_empty_dicts(v) for v in data]
    else:
        return data


class BaseModel(PydanticBaseModel):
    def __init__(__pydantic_self__, **data: Any):
        cleaned = remove_empty_dicts(data)
        super().__init__(**cleaned)


class TemplateMixin:
    id: Optional[str] = Field(default=None, description="ID in the graph")
    label: Optional[str] = Field(default=None, description="Label in the graph")

    @classmethod
    def template(cls, pid: bool = False) -> dict:
        def unwrap_type(typ):
            origin = get_origin(typ)
            if origin is Union:
                return next(
                    (arg for arg in get_args(typ) if arg is not type(None)), None
                )
            return typ

        template = {}
        for name, field in cls.model_fields.items():
            if not pid and name in {"pid", "label"}:
                continue

            typ = unwrap_type(field.annotation)

            if isinstance(typ, type) and issubclass(typ, BaseModel):
                if hasattr(typ, "template"):
                    template[name] = typ.template(pid=pid)
                else:
                    template[name] = {}

            elif get_origin(typ) is DataProperty:
                template[name] = {"value": None, "pid": None, "unit": None}

            elif field.default is not None:
                template[name] = field.default

            elif callable(field.default_factory):
                template[name] = field.default_factory()

            else:
                template[name] = None

        return template


class DataProperty(BaseModel, Generic[T], TemplateMixin):
    value: Optional[T] = None
    pid: Optional[str] = Field(default=None, description="PID")
    unit: Optional[str] = Field(default=None, description="Unit of measure")
    id: Optional[str] = Field(default=None, description="ID in the graph")
    label: Optional[str] = Field(default=None, description="Label in the graph")
    basename: Optional[str] = Field(
        default=None, description="Basename for file-based data properties"
    )


class Activity(BaseModel, TemplateMixin):
    pid: Optional[str] = Field(default=None, description="PID of the activity")

    initial_sample: Optional[str] = Field(
        default=None, description="ID of the initial sample in the graph"
    )
    final_sample: Optional[str] = Field(
        default=None, description="ID of the final sample in the graph"
    )


class RDFMixin:
    def to_graph(self, graph, sample):
        raise NotImplementedError

    @classmethod
    def from_graph(cls, graph, sample):
        raise NotImplementedError
