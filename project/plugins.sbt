resolvers += Resolver.url("fix-sbt-plugin-releases", url("http://dl.bintray.com/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)

addSbtPlugin("com.typesafe.sbt"  % "sbt-git"       % "0.9.3")
addSbtPlugin("com.github.gseitz" % "sbt-release"   % "1.0.6")
addSbtPlugin("com.eed3si9n"      % "sbt-assembly"  % "0.14.6")
addSbtPlugin("org.scoverage"     % "sbt-scoverage" % "1.5.1")
addSbtPlugin("org.xerial.sbt"    % "sbt-sonatype"  % "2.2")
addSbtPlugin("com.jsuereth"      % "sbt-pgp"       % "1.1.0")
addSbtPlugin("net.virtual-void"  % "sbt-dependency-graph" % "0.9.0")
addSbtPlugin("com.thoughtworks.sbt-api-mappings" % "sbt-api-mappings" % "2.0.0")

