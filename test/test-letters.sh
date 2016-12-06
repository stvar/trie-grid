#!/bin/bash

# Copyright (C) 2016  Stefan Vargyas
# 
# This file is part of Trie-Grid.
# 
# Trie-Grid is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Trie-Grid is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Trie-Grid.  If not, see <http://www.gnu.org/licenses/>.

#
# File generated by a command like:
# $ gen-test -T letters
#

[[ "$1" =~ ^-u[0-9]+$ ]] &&
u="${1:2}" ||
u=""

diff -u$u -L letters.old <(echo \
'$ test -x ./trie
$ print() { printf '\''%s\n'\'' "$@"; }
$ letters() { python -c $'\''\nfrom sys import stdout\nfor i in range(ord("a"),ord("z")+1):\n\tfor j in range(ord("a"),i+1):\n\t\tstdout.write(chr(j))\n\tstdout.write("\\n")\n'\''; }
$ letters|sort -c
$ letters|./trie -u -w -D
{
.   '\''a'\'': "a" {
.   .   '\''b'\'': "ab" {
.   .   .   '\''c'\'': "abc" {
.   .   .   .   '\''d'\'': "abcd" {
.   .   .   .   .   '\''e'\'': "abcde" {
.   .   .   .   .   .   '\''f'\'': "abcdef" {
.   .   .   .   .   .   .   '\''g'\'': "abcdefg" {
.   .   .   .   .   .   .   .   '\''h'\'': "abcdefgh" {
.   .   .   .   .   .   .   .   .   '\''i'\'': "abcdefghi" {
.   .   .   .   .   .   .   .   .   .   '\''j'\'': "abcdefghij" {
.   .   .   .   .   .   .   .   .   .   .   '\''k'\'': "abcdefghijk" {
.   .   .   .   .   .   .   .   .   .   .   .   '\''l'\'': "abcdefghijkl" {
.   .   .   .   .   .   .   .   .   .   .   .   .   '\''m'\'': "abcdefghijklm" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''n'\'': "abcdefghijklmn" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''o'\'': "abcdefghijklmno" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''p'\'': "abcdefghijklmnop" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''q'\'': "abcdefghijklmnopq" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''r'\'': "abcdefghijklmnopqr" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''s'\'': "abcdefghijklmnopqrs" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''t'\'': "abcdefghijklmnopqrst" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''u'\'': "abcdefghijklmnopqrstu" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''v'\'': "abcdefghijklmnopqrstuv" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''w'\'': "abcdefghijklmnopqrstuvw" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''x'\'': "abcdefghijklmnopqrstuvwx" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''y'\'': "abcdefghijklmnopqrstuvwxy" {
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   '\''z'\'': "abcdefghijklmnopqrstuvwxyz"
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   .   }
.   .   .   .   .   .   .   }
.   .   .   .   .   .   }
.   .   .   .   .   }
.   .   .   .   }
.   .   .   }
.   .   }
.   }
}
$ diff -u -Lcompact <(letters|./trie -u -c) -Lwide <(letters|./trie -u -w|ssed -R '\''s/\x27/\x22/g'\'') >/dev/null'
) -L letters.new <(
echo '$ test -x ./trie'
test -x ./trie 2>&1 ||
echo 'command failed: test -x ./trie'

echo '$ print() { printf '\''%s\n'\'' "$@"; }'
print() { printf '%s\n' "$@"; } 2>&1 ||
echo 'command failed: print() { printf '\''%s\n'\'' "$@"; }'

echo '$ letters() { python -c $'\''\nfrom sys import stdout\nfor i in range(ord("a"),ord("z")+1):\n\tfor j in range(ord("a"),i+1):\n\t\tstdout.write(chr(j))\n\tstdout.write("\\n")\n'\''; }'
letters() { python -c $'\nfrom sys import stdout\nfor i in range(ord("a"),ord("z")+1):\n\tfor j in range(ord("a"),i+1):\n\t\tstdout.write(chr(j))\n\tstdout.write("\\n")\n'; } 2>&1 ||
echo 'command failed: letters() { python -c $'\''\nfrom sys import stdout\nfor i in range(ord("a"),ord("z")+1):\n\tfor j in range(ord("a"),i+1):\n\t\tstdout.write(chr(j))\n\tstdout.write("\\n")\n'\''; }'

echo '$ letters|sort -c'
letters|sort -c 2>&1 ||
echo 'command failed: letters|sort -c'

echo '$ letters|./trie -u -w -D'
letters|./trie -u -w -D 2>&1 ||
echo 'command failed: letters|./trie -u -w -D'

echo '$ diff -u -Lcompact <(letters|./trie -u -c) -Lwide <(letters|./trie -u -w|ssed -R '\''s/\x27/\x22/g'\'') >/dev/null'
diff -u -Lcompact <(letters|./trie -u -c) -Lwide <(letters|./trie -u -w|ssed -R 's/\x27/\x22/g') >/dev/null 2>&1 ||
echo 'command failed: diff -u -Lcompact <(letters|./trie -u -c) -Lwide <(letters|./trie -u -w|ssed -R '\''s/\x27/\x22/g'\'') >/dev/null'
)

