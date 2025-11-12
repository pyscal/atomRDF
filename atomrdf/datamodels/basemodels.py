from typing import List, Optional, Union, Generic, TypeVar, get_args, get_origin, Any
from pydantic import BaseModel as PydanticBaseModel
from pydantic import Field

T = TypeVar("T")


def remove_empty_dicts(data):
    if isinstance(data, dict):
        cleaned = {}
        for k, v in data.items():
            cleaned_v = remove_empty_dicts(v)
            # Skip empty dicts
            if isinstance(cleaned_v, dict) and cleaned_v == {}:
                continue
            # Skip None values
            if cleaned_v is None:
                continue
            # Skip empty lists
            if isinstance(cleaned_v, list) and cleaned_v == []:
                continue
            # Keep everything else
            cleaned[k] = cleaned_v
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
            # Skip internal fields
            if name in {"graph", "id"}:
                continue
            # Skip pid unless explicitly requested
            if not pid and name == "pid":
                continue

            typ = unwrap_type(field.annotation)

            # Handle List types first
            if get_origin(typ) is list:
                args = get_args(typ)
                if args:
                    # Get the inner type (e.g., SomeModel from List[SomeModel])
                    inner_typ = unwrap_type(args[0])
                    # Check against PydanticBaseModel to catch all BaseModel subclasses
                    if isinstance(inner_typ, type) and issubclass(
                        inner_typ, PydanticBaseModel
                    ):
                        # Create a list with one template example
                        if hasattr(inner_typ, "template"):
                            template[name] = [inner_typ.template(pid=pid)]
                        else:
                            template[name] = [{}]
                    else:
                        # For primitive types, just empty list
                        template[name] = []
                else:
                    template[name] = []

            # Handle BaseModel types (including those with default=None)
            # Check against PydanticBaseModel to catch both custom BaseModel and direct Pydantic subclasses
            elif isinstance(typ, type) and issubclass(typ, PydanticBaseModel):
                if hasattr(typ, "template"):
                    template[name] = typ.template(pid=pid)
                else:
                    template[name] = {}

            # Handle DataProperty types
            elif get_origin(typ) is DataProperty:
                # DataProperty has template() method, use it
                if hasattr(typ, "template"):
                    template[name] = typ.template(pid=pid)
                else:
                    # Fallback if for some reason template is not available
                    template[name] = {"value": None, "pid": None, "unit": None}

            # Handle fields with non-None defaults (but only for non-BaseModel types)
            elif field.default is not None:
                template[name] = field.default

            # Handle default_factory
            elif (
                hasattr(field, "default_factory")
                and field.default_factory is not None
                and callable(field.default_factory)
            ):
                template[name] = field.default_factory()

            # Everything else gets None
            else:
                template[name] = None

        return template

    @classmethod
    def save_template(
        cls, filepath: str, format: str = "json", pid: bool = False, indent: int = 2
    ) -> None:
        """
        Save the template to a file in JSON or YAML format.

        Parameters
        ----------
        filepath : str
            Path where the template file will be saved
        format : str, optional
            Output format, either 'json' or 'yaml' (default: 'json')
        pid : bool, optional
            Whether to include PID fields in the template (default: False)
        indent : int, optional
            Number of spaces for JSON indentation (default: 2, only used for JSON format)

        Examples
        --------
        >>> AtomicScaleSample.save_template('template.json', format='json')
        >>> AtomicScaleSample.save_template('template.yaml', format='yaml')
        >>> AtomicScaleSample.save_template('template.json')  # defaults to JSON
        """
        template = cls.template(pid=pid)

        if format.lower() == "json":
            import json

            with open(filepath, "w") as f:
                json.dump(template, f, indent=indent)
        elif format.lower() == "yaml":
            import yaml

            with open(filepath, "w") as f:
                yaml.dump(template, f, default_flow_style=False, sort_keys=False)
        else:
            raise ValueError(f"Unsupported format: '{format}'. Use 'json' or 'yaml'.")


class DataProperty(BaseModel, Generic[T], TemplateMixin):
    value: Optional[T] = None
    pid: Optional[str] = Field(default=None, description="PID")
    unit: Optional[str] = Field(default=None, description="Unit of measure")
    basename: Optional[str] = Field(
        default=None, description="Basename for file-based data properties"
    )


class RDFMixin:
    def to_graph(self, graph, sample):
        raise NotImplementedError

    @classmethod
    def from_graph(cls, graph, sample):
        raise NotImplementedError
