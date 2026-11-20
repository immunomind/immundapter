

# Set the default output directory for file-based runs

[**Source code**](https://github.com/immunomind/immundapter/tree/main/R/output-dir.R#L11)

## Description

Save a directory for file-based runs across R sessions. An explicit
<code>output_dir</code> passed to
<code style="white-space: pre;">$run()</code> takes precedence. The
directory is created when a run first uses it.

## Usage

<pre><code class='language-R'>set_output_dir(path = NULL)
</code></pre>

## Arguments

<table role="presentation">
<tr>
<td style="white-space: nowrap; font-family: monospace; vertical-align: top">
<code id="path">path</code>
</td>
<td>
Directory to save. In an interactive session, omit <code>path</code> to
enter and confirm a directory at the console.
</td>
</tr>
</table>

## Value

Invisibly, the absolute path that was saved, or <code>NULL</code> if
cancelled.
