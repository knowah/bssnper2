#/bin/sh
HTSLIB=../../software/htslib
gcc -g -Wall -o bssnper2 main.c bssnper2.c genotyping.c uncalled_pos.c -I"$HTSLIB" -L"$HTSLIB" -lhts -lm -lz -lpthread -lbz2 -llzma -lcurl -lcrypto