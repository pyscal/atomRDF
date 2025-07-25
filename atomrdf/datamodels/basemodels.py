from typing import List, Optional, Union, Generic, TypeVar, get_args, get_origin
from pydantic import BaseModel, Field

T = TypeVar("T")


class TemplateMixin:
    @classmethod
    def template(cls) -> dict:
        def unwrap_type(typ):
            origin = get_origin(typ)
            if origin is Union:
                return next(
                    (arg for arg in get_args(typ) if arg is not type(None)), None
                )
            return typ

        template = {}
        for name, field in cls.model_fields.items():
            typ = unwrap_type(field.annotation)

            if isinstance(typ, type) and issubclass(typ, BaseModel):
                template[name] = typ.template() if hasattr(typ, "template") else {}

            elif get_origin(typ) is DataProperty:
                template[name] = {"value": None, "pid": None, "unit": None}
            else:
                template[name] = None

        return template


class DataProperty(BaseModel, Generic[T]):
    value: Optional[T] = None
    pid: Optional[str] = Field(default=None, description="PID")
    unit: Optional[str] = Field(default=None, description="Unit of measure")
