#
# Makefile for the third program
#
CP=jars/Graph.jar:jars/SteinerTreeTester.jar 
RUN_CP=SteinerTreeTester.jar:student.jar:Graph.jar
build:
	javac -cp $(CP) -d bin src/student/SteinerTree.java
	cd bin; \
	jar cf ../jars/student.jar student
run:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester
test:
