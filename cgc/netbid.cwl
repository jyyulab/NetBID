{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "id": "cjrash/netbid/netbid/11",
    "baseCommand": [],
    "inputs": [
        {
            "id": "rscript",
            "type": "File",
            "inputBinding": {
                "shellQuote": false,
                "position": 0
            }
        },
        {
            "id": "expression_set",
            "type": "File",
            "inputBinding": {
                "prefix": "-e",
                "shellQuote": false,
                "position": 1
            },
            "label": "Gene Expression Matrix"
        },
        {
            "id": "tf",
            "type": "File",
            "inputBinding": {
                "prefix": "-t",
                "shellQuote": false,
                "position": 1
            },
            "label": "Transcription Factor Network"
        },
        {
            "id": "sig",
            "type": "File",
            "inputBinding": {
                "prefix": "-s",
                "shellQuote": false,
                "position": 1
            },
            "label": "Signaling Network"
        },
        {
            "id": "metadata",
            "type": "File",
            "inputBinding": {
                "prefix": "-m",
                "shellQuote": false,
                "position": 3
            },
            "label": "Sample Grouping File"
        },
        {
            "sbg:toolDefaultValue": "project",
            "id": "project_name",
            "type": "string?",
            "inputBinding": {
                "prefix": "-p",
                "shellQuote": false,
                "position": 4
            }
        }
    ],
    "outputs": [
        {
            "id": "output",
            "type": "Directory?",
            "outputBinding": {
                "glob": "NetBID_*",
                "loadListing": "deep_listing"
            }
        },
        {
            "id": "netbid_log",
            "type": "stdout",
            "outputBinding": {
                "glob": "$(inputs.expression_set.nameroot).log"
            }
        }
    ],
    "label": "netbid",
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "DockerRequirement",
            "dockerPull": "cgc-images.sbgenomics.com/cjrash/stjude/netbid:latest"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "stdout": "$(inputs.expression_set.nameroot).log",
    "sbg:projectName": "netbid",
}
