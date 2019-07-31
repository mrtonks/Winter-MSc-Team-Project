macOS and Linux Systems
=======
bash test.sh

it can install docker and run it automatically


Windows Systems
=======
bash test_windows.sh

IMPORTANT NOTE: Before executing the command:

1. Make sure the project folder is on a location within "C:\Users\"
2. Locate the folder inside the project "docker\Pipeline" 
3. Locate and open the file inside the project folder "docker\test_windows.sh" on a text editor (i.e. Notepad)
4. Change the line 5 with the physical path of the "Pipeline" folder
	Change -> /c/Users/inspiron/Pipeline
	By -> /c/Users/<project_folder>/docker/Pipeline
	
	Line 5 should look like: docker run --rm -it -v /c/Users/<project_folder>/docker/Pipeline:/home/pipeline/Pipeline \
5. Save the file
6. Run the command: bash test_windows.sh