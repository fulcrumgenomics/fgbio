resolvers += Resolver.url("fix-sbt-plugin-releases", url("https://dl.bintray.com/sbt/sbt-plugin-releases"))(Resolver.ivyStylePatterns)

addDependencyTreePlugin

addSbtPlugin("com.typesafe.sbt"  % "sbt-git"       % "1.0.0")
addSbtPlugin("com.github.gseitz" % "sbt-release"   % "1.0.10")
addSbtPlugin("com.eed3si9n"      % "sbt-assembly"  % "1.2.0")
addSbtPlugin("org.scoverage"     % "sbt-scoverage" % "1.6.1")
addSbtPlugin("com.github.sbt"    % "sbt-pgp"       % "2.3.1")
addSbtPlugin("com.thoughtworks.sbt-api-mappings" % "sbt-api-mappings" % "3.0.0")
