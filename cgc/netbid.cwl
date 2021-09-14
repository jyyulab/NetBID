{
    "class": "CommandLineTool",
    "cwlVersion": "v1.2",
    "baseCommand": [],
    "inputs": [
        {
            "id": "expression_set",
            "type": "File",
            "inputBinding": {
                "prefix": "-e",
                "shellQuote": false,
                "position": 1
            },
            "label": "Gene Expression Matrix",
            "doc": "comma-delimited expression matrix file with columns as samples, rows as genes."

        },
        {
            "id": "tf",
            "type": "File",
            "inputBinding": {
                "prefix": "-t",
                "shellQuote": false,
                "position": 1
            },
            "label": "Transcription Factor Network",
            "doc": "file with each row an edge from the TF network constructed using SJARACNe (https://github.com/jyyulab/SJARACNe)"

        },
        {
            "id": "sig",
            "type": "File",
            "inputBinding": {
                "prefix": "-s",
                "shellQuote": false,
                "position": 1
            },
            "label": "Signaling Network",
            "doc": "file with each row an edge from the SIG network constructed using SJARACNe (https://github.com/jyyulab/SJARACNe)"

        },
        {
            "id": "metadata",
            "type": "File",
            "inputBinding": {
                "prefix": "-m",
                "shellQuote": false,
                "position": 3
            },
            "label": "Sample Grouping File",
            "doc": "comma-delimited file with two columns: sample and group."
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
    "doc": "# Description\n\nNetBID is a data-driven system biology pipeline using a data-driven network-based Bayesian inference approach to find drivers from transcriptomics, proteomics, or phosphoproteomics data. The drivers can be either transcription factors (TF) or signaling factors (SIG).\n\nNetBID2 has the following key steps to perform hidden driver analysis:\n1.\tActivity calculation of drivers based on drivers’ regulons from a pre-built or user-provided SJARACNe network;\n2.\tDiscovery of differential expressed genes and differential activated drivers;\n3.\tGeneration of the master table for drivers;\n4.\tVisualizing drivers with significance profiles and target genes.\n\n# Inputs and outputs of NetBID workflow\n## Inputs:\n*\tExpression matrix - comma-delimited expression matrix file with columns as samples, rows as genes.\n*\tMetadata file - comma-delimited file with two columns: sample and group.\n*\tSignaling (SIG) network - file with each row an edge from the SIG network constructed using SJARACNe (https://github.com/jyyulab/SJARACNe)\n*\tTranscription factor (TF) network - file with each row an edge from the TF network constructed using SJARACNe (https://github.com/jyyulab/SJARACNe)\n\n## Outputs:\n*\tExcel file with differential expressed genes and differential activated drivers\n*\tPicture file visualizing drivers with significance profiles\n\n# Common issues\n*\tThe first row and the first column of the expression matrix file must be sample names and gene names, respectively.\n*\tThe metadata file must have at least two sample groups in the 2nd column.",

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
            "dockerPull": "cgc-images.sbgenomics.com/stjude/netbid:latest"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "stdout": "$(inputs.expression_set.nameroot).log",
    "sbg:projectName": "netbid",
}
