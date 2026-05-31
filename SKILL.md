# Pipeline Skill: Task Logging

Use this skill at the end of any task to produce a consistent 
condensed summary log.

## When to use
Append "Follow SKILL.md logging format." 
to any prompt to activate this skill.

## Log format
Write a plain text file with these sections, in order:

TASK
One sentence describing what was asked.

FOUND
What was discovered — relevant code, data, or state. Be specific: 
file name, function name, line number if applicable.

CHANGED
Exact lines, functions, or files modified. If nothing changed, say so.

RESULT
Confirmation that the expected outcome was achieved, or description 
of what blocked it.

COMMIT
Hash and full commit message. If no commit was made, say so.

## Rules
- Maximum 20 lines total
- No markdown formatting inside the log file — plain text only
- Save in project root as session.log — always this name, overwritten each session
- Never summarise what you plan to do — only what was done
