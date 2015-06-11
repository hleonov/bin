#!/bin/csh -f

evolution --force-shutdown
find ~/.evolution/mail/local -name "*ev-summary" -exec rm -f {} \;
