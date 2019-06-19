build-java:
	@cd hivfacts-java; ./gradlew build

release-java:
	@cd hivfacts-java; ./gradlew build bintrayUpload

build-python:
	@cd hivfacts-python; ./setup.py copy_data sdist bdist_wheel

release-python:
	@cd hivfacts-python; ./setup.py copy_data sdist bdist_wheel upload

.PHONY: release-java
