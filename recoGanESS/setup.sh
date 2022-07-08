#!/usrbin/env bash

export RGANESS=$PWD
export GRESDIR=$PWD/gres
export PATH=$RGANESS/bin:$PATH
export PATH=$(tr ":" "\n" <<<"$PATH" | grep -Fxv "$ICTDIR/bin" | paste -sd:)
export PYTHONPATH=$RGANESS:$PYTHONPATH

