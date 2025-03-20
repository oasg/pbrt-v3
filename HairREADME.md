# Hair Model

## Adding a New Model
Hair models are located in ```src/materials/mhair_new.h```, mainly inheriting from the hairbsdf class. To modify the reflection model, you only need to modify the ```Spectrum f(const Vector3f &wo, const Vector3f &wi) const;``` function.

To add a new model, you need to modify ```src/core/api.cpp``` to add a new material, include the header file, and add a fragment for pbrt model parsing, so that the compiler knows the definition of the new model.

The ```class SingBrdf``` class is a singleton that reads the table, ensuring that it is only loaded once when the program runs.

The multiple reflection part of the original hair model is currently retained. If you don't need to consider scattering effects, you can comment out this part, or directly modify the parameters of the scene.
## Scene File
Add a new field in the scene file to receive the location of the table. You can modify it here when changing the model.
```
"string table_path" ["../table/HairMultilayerDamagedTiltPerlinModel/incidenceLayer/(10nm,1200cell)/TM/Ns/Reflection/"]
```

## Other Tools
In the ```tools``` directory, using some of the built-in features of pbrt, a color transformation program is written. Note that its function is to visualize the spectral data in the table into an image.

In the ```table/visualize``` directory, there are some scripts to visualize the table.