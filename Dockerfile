FROM eclipse-temurin:11
RUN adduser --system --group --shell /bin/sh dockeruser \
 && mkdir /home/dockeruser/bin
RUN mkdir /opt/app
COPY target/biozernike-1.0.0-SNAPSHOT-jar-with-dependencies.jar /opt/app/ligzernike.jar
USER dockeruser
ENTRYPOINT ["java", "-jar", "/opt/app/ligzernike.jar"]

