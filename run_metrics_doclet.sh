#!/bin/bash

###############################################################################
# Simple shell script to generate metrics documentation with scaladoc
###############################################################################

# JAVA_OPTS="-agentlib:jdwp=transport=dt_socket,server=y,suspend=n,address=7777"

sources=$(find src/main/scala -name \*.scala)
cp=$(sbt -Dsbt.log.noformat=true "export runtime:fullClasspath" 2> /dev/null | fgrep -v '[info]')

scaladoc -toolcp $cp -d target -doc-generator com.fulcrumgenomics.internal.FgMetricsDoclet $sources