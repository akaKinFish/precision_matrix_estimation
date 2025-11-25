---
type: "manual"
---

# Coding and Debugging Workflow
## patience
there are no time limits, no token limits, some times after some short way with scripts, we need check/fix the result manualy
# document
-- don't generate too much md file, too redundant
-- don't stop and wait user's feedback before tasks finish
everything in English

# Don't write any md/txt or any other files for summary
- Don't write md/txt or any other files for summary, it is waste of token, since you never use them again.
# better not use unicode or emoji in the code, aovid error
# CODE STYLE
the code should High Cohesion, Low Coupling, Modularity

"TODO.md" is the task tracker, please read and update every chat. It is very big, so Operate it parts by parts?

"doc\janlag_eeg_application_plan.md" is the overall project detail.

please keep all important folder clean.


important docs (paper pdf or ref code, plan md et. al.)  in ./doc/
test code please put in ./dev/
please when you generate summary put ./dev/doc/
output result put ./result/
output figure put ./fig/
useful data put ./data/
temp file put ./temp/
no longer useful and wrong things put ./backup/

3rd part package in ./externel/ don't add all to path, add useful folder
make the folder clean, move old files to ./dev/old/
 

don't keep files or create files like "_v1", "_v2", "_fixed" et. al. move the old to ./backup/ and then only keep one useful, keep folder clearn


FOR MATLAB CODE: PLEASE USE MATLAB VECTORIZATION TO ACCELERATE

## General Principles
- Treat every task as two deliverables: the implementation and the verified evidence (logs, screenshots, test results) that it works.
- Prefer reproducible command-line steps (PowerShell or CMD) for building, testing, linting, and inspecting artifacts.
- Record any environment setup or configuration changes so the user can rerun the same commands later.
- Use MCP tools whenever they help gather documentation, explore resources, or automate browser interactions; note which tool you invoked and why.

## Unit Testing
- Add necessary unit test code for each function.
- Tests should not only check the type and size of data but also verify the value of critical data against theoretical expectations.

## MCP Usage
- Available servers: `context7`, `everything`, `fetch`, `filesystem`, `mcp-deepwiki`, `memory`, `open-websearch`, `playwright`, `sequential-thinking`, `time`.
- Typical commands:
  - Context7 docs: `mcp__mcp-router__get-library-docs` with a library ID (e.g., `/vercel/next.js`).
  - Static resources: `read_mcp_resource` for `resource://everything/test://static/resource/...`.
  - Web fetch: `mcp__mcp-router__fetch <url>` for raw HTTP content.
  - Allowed paths and metadata: `mcp__mcp-router__list_allowed_directories`, `mcp__mcp-router__get_file_info`.
  - Web search: `mcp__mcp-router__search`.
  - Browser automation: `mcp__mcp-router__browser_*` (navigate, snapshot, interactions).
  - Knowledge graph memory: `mcp__mcp-router__read_graph` (and corresponding write/update calls when needed).
  - Time utilities: `mcp__mcp-router__get_current_time`.
  - DeepWiki content: `mcp__mcp-router__deepwiki_fetch`.
  - Structured reasoning: `mcp__mcp-router__sequentialthinking`.
- Mention MCP usage in your response when it materially contributed to the result so the user can trace the workflow.
## Media data extraction
Use 'external\gemini-api\process_pdf.py' to use gemini-api to process PDFs or other image, audio files. it is multimodal LLM, embedded all medias in the same common latent space, which understands the context of the image, audio, and can answer questions about it perfectly.

## Debugging and Verification Checklist
- After writing code, run the relevant executable or test suite; capture console output or generated artifacts.
- Inspect error logs and warnings even if the command succeeds; document any non-blocking issues.
- For scripts producing files, confirm their presence, size, and key contents when applicable.
- Summarize verification steps in the final update so the user knows exactly how you confirmed correctness.

### Visual Verification
- For visual checking later by the programmer/user, we can leave some plot of essential variables and activate these by some flags.

### Debugging and Iteration
- If the debug/verification has a problem, please fix the problem until you go to the next step. It must be fixed.

## MATLAB Workflow
1. Write the required `.m` files with clear structure; include `saveas` for figures.
2. Execute via `matlab -batch "run('script_name.m')"` from the working directory.
3. Collect console output; if files were produced (`.mat`, `.csv`, images), verify them with follow-up commands (`dir`, `type`, etc.).
4. Report success or highlight issues; attach figure filenames or data summaries.

### Example MATLAB Workflow
- **Request:** Create a sine-wave plot.
  1. Author `plot_sine.m` producing `sine_wave.png`.
  2. Run `matlab -batch "run('plot_sine.m')"`。
  3. Confirm `sine_wave.png` exists (e.g., `dir sine_wave.png`) and note key details.
  4. Inform the user that execution succeeded and the image was generated.
  5. MATLAB private folders: Functions in utility/private/ can ONLY be called by MATLAB functions (not scripts) in the parent folder utility/ - this is a MATLAB mechanism, not a naming convention.

## Python Workflow
1. Check for `.venv`; if absent, run `python -m venv .venv`.
2. Activate the environment:
   - CMD: `.\.venv\Scripts\activate`
   - PowerShell: `.\.venv\Scripts\Activate.ps1`
   - Linux/macOS: `source .venv/bin/activate`
3. Use `uv` for package management: `uv pip install <package>` or `uv pip install -r requirements.txt`.
4. Implement the Python code inside the project workspace.
5. Execute scripts/tests with the venv interpreter (e.g., `python script.py`, `python -m pytest`).
6. Capture outputs, logs, and any generated files; note imported packages that confirm dependency installation.

### Example Python Workflow
- **Request:** Script using `numpy`.
  1. Ensure `.venv` exists; create if needed.
  2. Activate the venv.
  3. Install `numpy` via `uv pip install numpy`.
  4. Write the script (e.g., `compute_stats.py`).
  5. Run `python compute_stats.py`; record the printed statistics.
  6. Share the results and confirm `numpy` executed properly.

## R Workflow
1. Place `.R` scripts in the workspace and ensure paths are correct.
2. Run scripts using the fixed executable: `"D:\Software\R\R-4.5.1\bin\R.exe" CMD BATCH script_name.R`.
3. Inspect the generated `.Rout` file for output and errors (`type script_name.Rout` or equivalent).
4. Report computed values and any warnings encountered.

### Example R Workflow
- **Request:** Mean calculation.
  1. Author `calculate_mean.R` printing the mean.
  2. Execute with `"D:\Software\R\R-4.5.1\bin\R.exe" CMD BATCH calculate_mean.R`.
  3. Read `calculate_mean.Rout`, extract the printed mean, and return it with confirmation of success.

## Closing Your Response
- Summarize the code changes or scripts created, referencing file paths.
- List the exact commands run for verification and highlight their outcomes.
- Flag any follow-up tasks the user may want to run (e.g., longer tests, deployment steps) and mention if any verification was skipped due to environment limits.

