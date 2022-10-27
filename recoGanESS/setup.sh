#!/usrbin/env bash

export RGANESS=$PWD
export GRESDIR=$PWD/gres
export PATH=$RGANESS/bin:$PATH
export AUXPATH=$(tr ":" "\n" <<<"$PATH" | grep -Fxv "$ICTDIR/bin")
export PATH=$(tr "\n" ":" <<< "$AUXPATH")
export PYTHONPATH=$RGANESS:$PYTHONPATH

