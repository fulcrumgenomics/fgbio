resolvers += Resolver.url("fix-sbt-plugin-releases", url("http://dl.bintray.com/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)

addSbtPlugin("com.typesafe.sbt"  % "sbt-git"       % "1.0.0")
addSbtPlugin("com.github.gseitz" % "sbt-release"   % "1.0.10")
addSbtPlugin("com.eed3si9n"      % "sbt-assembly"  % "0.15.0")
addSbtPlugin("org.scoverage"     % "sbt-scoverage" % "1.6.1")
addSbtPlugin("org.xerial.sbt"    % "sbt-sonatype"  % "3.9.4")
addSbtPlugin("com.jsuereth"      % "sbt-pgp"       % "2.0.1")
addSbtPlugin("net.virtual-void"  % "sbt-dependency-graph" % "0.9.2")
addSbtPlugin("com.thoughtworks.sbt-api-mappings" % "sbt-api-mappings" % "2.0.0")

