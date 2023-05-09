#!/usr/bin/env python3

import json
import os
import subprocess
import shlex

script = os.path.realpath(__file__)
script_dir = os.path.dirname(script)
workflows = json.load(open(os.path.join(script_dir, "workflows.json")))
wf_id_file = os.path.join(script_dir, "workflow-ids.json")

patch_cmd = "patch {input_wdl} {diff_file} -o {out_wdl}"
create_cmd = "aws omics create-workflow --name {wf_name} --description {wf_description} --definition-zip fileb://{wf_zip} --parameter-template file://{template} --engine WDL"

def main():
    os.chdir(os.path.join(script_dir, "..", "pipelines"))
    completed_wf_ids = {}

    for wf_name, wf in workflows.items():
        wdl = wf["wdl"]
        if "diff" in wf and wf["diff"]:  # Patch the wdl
            out_wdl = wf["diff"].replace(".diff", ".wdl")
            cmd = shlex.split(patch_cmd.format(input_wdl=wdl, diff_file=wf["diff"], out_wdl=out_wdl))
            subprocess.run(cmd, check=True)
            wdl = out_wdl

        # Zip the workflow
        wf_zip = wdl.replace(".wdl", ".zip")
        subprocess.run(shlex.split(f"zip {wf_zip} {wdl}"), check=True)

        # Create the workflow
        cmd = shlex.split(create_cmd.format(
            wf_name=wf["name"],
            wf_description='"' + wf["description"] + '"',
            wf_zip=wf_zip,
            template=wf["template"],
        ))
        p = subprocess.run(cmd, check=True, text=True, stdout=subprocess.PIPE)
        p_json = json.loads(p.stdout)

        completed_wf_ids[wf_name] = p_json["id"]

        os.unlink(wf_zip)
        if "diff" in wf and wf["diff"]:
            os.unlink(wdl)

    with open(wf_id_file, 'w') as fh:
        json.dump(completed_wf_ids, fh)

if __name__ == "__main__":
    main()
