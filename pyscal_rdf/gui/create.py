import ipywidgets as widgets
from ipywidgets import Layout
from pyscal_rdf.gui.themes import themes
import panel as pn
pn.extension()

def output(theme="teal"):
    """
    Create an ouput widget
    """
    output = widgets.Output(layout={'border': '1px solid %s'%themes[theme]["border"]})
    return output

def dropdown(description, options, value=None, disabled=False, tooltip="Select value", theme="teal"):
    if value is None:
        value = options[0]
        
    dropdown = widgets.Dropdown(
        options = options,
        value = value,
        description = description,
        description_long = description,
        disabled=disabled,
        tooltip=tooltip,
    )    
    return dropdown

def button(description, disabled=False, tooltip="Click me", theme="teal"):
    button = widgets.Button(
        description=description,
        disabled=disabled,
        button_style='',
        tooltip=tooltip,
    )
    button.style.button_color = themes[theme]["button"]
    return button

def checkbox(description, value=True, theme="teal"):
    checkbox = widgets.Checkbox(
        value=value,
        description=description,
        disabled=False
    )
    return checkbox 

def textbox(description, value, dtype, theme="teal", disabled=False, tooltip="Enter value"):
    if dtype == "int":
        inttext = widgets.IntText(
            value=value,
            description=description,
            disabled=disabled,
            tooltip=tooltip,
        )
        return inttext
    elif dtype == "float":
        inttext = widgets.FloatText(
            value=value,
            description=description,
            disabled=disabled,
            tooltip=tooltip,
        )
        return inttext
    elif dtype == "text":
        inttext = widgets.Text(
            value=value,
            description=description,
            disabled=disabled,
            tooltip=tooltip,
        )
        return inttext 
    elif dtype == "textarea":
        inttext = widgets.Textarea(
            value=value,
            description=description,
            disabled=disabled,
            tooltip=tooltip,
            layout=Layout(width="auto", height="100%")
        )
        return inttext 

def header(text, tooltip="", theme="teal"):
    color = themes[theme]["header"]
    header = widgets.HTML(value = f"<b><font color='{color}'>{text}</b>",tooltip=tooltip)
    return header

def text(text, tooltip="", theme="teal"):
    color = themes[theme]["text"]
    header = widgets.HTML(value = f"<font color='{color}'>{text}", tooltip=tooltip)
    return header

def upload(theme="teal", tooltip="Select file to upload"):
    upload = widgets.FileUpload(
        accept='',
        multiple=False,
        tooltip=tooltip,
    )
    return upload

def download(filename, theme="teal"):
    download = pn.widgets.FileDownload(file=filename, 
                            embed=True)
    return download
    