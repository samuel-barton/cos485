#
## Makefile for the third program
#
CP=jars/Graph.jar:jars/SteinerTreeTester.jar:bin/
RUN_CP=SteinerTreeTester.jar:student.jar:Graph.jar
build:
	javac -cp $(CP) -d bin src/student/WeightedVertex.java
	javac -cp $(CP) -d bin src/student/SteinerTree.java
	cd bin; \
	jar cf ../jars/student.jar student
run:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester 
test1:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test1.txt
test2:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test2.txt
test3:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test3.txt
test4:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test4.txt
test5:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test5.txt
test6:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test6.txt
test7:
	cd jars; \
	java -cp $(RUN_CP) steinerTree.SteinerTreeTester test7.txt
