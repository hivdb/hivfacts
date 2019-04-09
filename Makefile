build-java:
	@cd hiv-aapcnt-java; ./gradlew build

release-java:
	@cd hiv-aapcnt-java; ./gradlew build bintrayUpload

build-python:
	@cd hiv-aapcnt-python; ./setup.py copy_aapcnt sdist bdist_wheel

release-python:
	@cd hiv-aapcnt-python; ./setup.py copy_aapcnt sdist bdist_wheel upload

.PHONY: release-java
