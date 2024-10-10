@echo off


rem Root OSGEO4W home dir to the same directory this script exists in
call "%~dp0\bin\o4w_env.bat"

set SAGA=C:/OSGeo4W/apps\saga
set SAGA_MLB=C:/OSGeo4W/apps\saga\tools

PATH=C:\OSGeo4W\apps\qgis-ltr\bin;C:\OSGeo4W\apps\grass\grass78\lib;C:\OSGeo4W\apps\grass\grass78\bin;C:\OSGeo4W\apps\Qt5\bin;C:\OSGeo4W\apps\Python39\Scripts;C:\OSGeo4W\bin;C:\Windows\System32;C:\Windows;C:\Windows\System32\wbem;C:\OSGeo4W\apps\saga;C:\OSGeo4W\apps\saga\tools

cmd.exe