
The core action loop is the following:

```r
# idata: ImmunData with your data

# Create a config for a tool, TCRdist3 in the example
conf <- conf_tcrdist3(cpu = 4)

# Create the tool instance to run the tool multiple times with the same config if needed
tool <- tool_tcrdist3(conf)

# Run the tool on the `input` - either ImmunData, file paths or data frames, depending on the tool
run <- tool$run(idata)

# Get the run results once it is complete
results <- run$results()

# Annotate the original ImmunData with clusters or other info from the tool results
idata_new <- results$annotate(idata)
```

Shorter version:

```r
idata_new <- conf_tcrdist3(cpu = 4) |> tool_tcrdist3()$run(idata)$results$annotate(idata)
```
