release-java:
	@cd hiv-aapcnt-java; ./gradlew build bintrayUpload

.PHONY: release-java
