#!/bin/bash

# This clears all the generated code,
# and the code-generator make files
export COMMAND="rm -r build/* *_codegen autogen_code"
echo $COMMAND
$COMMAND
