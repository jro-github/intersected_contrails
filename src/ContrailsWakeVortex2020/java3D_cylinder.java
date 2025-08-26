package ContrailsWakeVortex2020;

import com.sun.j3d.utils.geometry.Cylinder;
import com.sun.j3d.utils.universe.SimpleUniverse;

import javax.media.j3d.*;
import javax.vecmath.*;

import javax.swing.*;

import java.io.IOException;


public class java3D_cylinder {
    public java3D_cylinder() throws IOException {
        // Create a 3D canvas
        Canvas3D canvas3D = new Canvas3D(SimpleUniverse.getPreferredConfiguration());
        canvas3D.setSize(1000, 1000);
        
        // Set up SimpleUniverse
        SimpleUniverse universe = new SimpleUniverse(canvas3D);
        universe.getViewingPlatform().setNominalViewingTransform();
        
        // Set different viewpoints
        setViewpoint(universe,new Vector3f(1.0f, 1.0f, 10.0f));

        // Create the scene
        BranchGroup scene = createTwoCylinders();

        // Create and add a background
        scene.addChild(createBackground());

        // Add the scene to Universe
        universe.addBranchGraph(scene);

        // Set up the window
        JFrame frame = new JFrame("Java 3D Elliptical Cylinders");
        frame.setSize(1000, 1000);
        frame.add(canvas3D);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setVisible(true);
    }
    
    private void setViewpoint(SimpleUniverse universe, Vector3f position) {
        TransformGroup viewTransformGroup = universe.getViewingPlatform().getViewPlatformTransform();

        Transform3D viewTransform = new Transform3D();
        viewTransform.setTranslation(position);

        // Create a downward tilt rotation
        Transform3D rotation = new Transform3D();
        
        viewTransform.mul(rotation); // Combine position and rotation

        viewTransformGroup.setTransform(viewTransform);
    }

    private BranchGroup createBackground() {
        // Create a branch group for the background
        BranchGroup bgNode = new BranchGroup();

        // Create a background and set its color
        Background background = new Background(new Color3f(0.9f, 0.9f, 0.9f));
        
        // Set the application bounds of the background
        BoundingSphere bounds = new BoundingSphere(new Point3d(0.0, 0.0, 0.0), 100.0);
        background.setApplicationBounds(bounds);

        // Add the background to the branch group
        bgNode.addChild(background);

        return bgNode;
    }
    

    private BranchGroup createTwoCylinders() {
        // Create a branch group
        BranchGroup objRoot = new BranchGroup();

        // Create transparent appearances
        Appearance app1 = new Appearance();
        Appearance app2 = new Appearance();

        // Set transparency 1
        TransparencyAttributes ta1 = new TransparencyAttributes();
        ta1.setTransparencyMode(TransparencyAttributes.SCREEN_DOOR);
        ta1.setTransparency(0.6f); // 0.0 completely opaque, 1.0 fully transparent
        app1.setTransparencyAttributes(ta1);
        
        // Set transparency 2
        TransparencyAttributes ta2 = new TransparencyAttributes();
        ta2.setTransparencyMode(TransparencyAttributes.SCREEN_DOOR);
        ta2.setTransparency(0.4f); // 0.0 completely opaque, 1.0 fully transparent
        app2.setTransparencyAttributes(ta2);

        // Set outline 1
        PolygonAttributes pa1 = new PolygonAttributes();
        pa1.setPolygonMode(PolygonAttributes.POLYGON_LINE);
        pa1.setCullFace(PolygonAttributes.CULL_NONE);
        app1.setPolygonAttributes(pa1);      
        
        // Set outline 2
        PolygonAttributes pa2 = new PolygonAttributes();
        pa2.setPolygonMode(PolygonAttributes.POLYGON_LINE);
        pa2.setCullFace(PolygonAttributes.CULL_NONE);
        app2.setPolygonAttributes(pa2);  
        
        // Set cylinder color to red 1
        Material material1 = new Material();
        material1.setDiffuseColor(new Color3f(1.0f, 0.0f, 0.0f));
        material1.setAmbientColor(new Color3f(1.0f, 0.0f, 0.0f));
        app1.setMaterial(material1);
        
        // Set cylinder color to blue 2
        Material material2 = new Material();
        material2.setDiffuseColor(new Color3f(0.0f, 0.0f, 1.0f));
        material2.setAmbientColor(new Color3f(0.0f, 0.0f, 1.0f));
        app2.setMaterial(material2);

        // Create and add cylinder 1
        TransformGroup tgCylinder1 = new TransformGroup();              
        // Create a transform object and set rotation 1
        Transform3D transform1 = new Transform3D();                
        transform1.rotZ(Math.PI / 2);  // Rotate 90 degrees around the Z axis
        transform1.setScale(new Vector3d(1.0, 3.0, 2.0)); 
        tgCylinder1.setTransform(transform1);
        
        // Create a cylinder
        Cylinder cylinder1 = new Cylinder(0.5f, 1.0f, app1);        
        tgCylinder1.addChild(cylinder1);
        
        // Create and add cylinder 2
        TransformGroup tgCylinder2 = new TransformGroup();              
        // Create a transform object and set rotation 2
        Transform3D transform2 = new Transform3D();     
        transform2.rotZ(Math.PI / 2);  // Rotate 90 degrees around the Z axis
        
        Transform3D rotateY = new Transform3D();
        rotateY.rotX(Math.PI / 2);  // Rotate 45 degrees around the Y axis
        transform2.mul(rotateY);
        transform2.setScale(new Vector3d(1.0, 5.0, 2.0)); 
        tgCylinder2.setTransform(transform2);
        
        // Create a cylinder
        Cylinder cylinder2 = new Cylinder(0.25f, 0.5f, app2);        
        tgCylinder2.addChild(cylinder2);
                
        // Add the cylinders 1 and 2 to the branch group
        objRoot.addChild(tgCylinder1);        
        objRoot.addChild(tgCylinder2);  
        
        // Add a coordinate system
        addCoordinateSystem(objRoot);

        return objRoot;
    }

    private void addCoordinateSystem(BranchGroup objRoot) {
        // Create an appearance for the coordinate axes
        //Appearance axisAppearance = new Appearance();
        //ColoringAttributes ca = new ColoringAttributes();
        //ca.setColor(1.0f, 0.0f, 0.0f); // X axis - Red
        //axisAppearance.setColoringAttributes(ca);
        //ca.setColor(0.0f, 1.0f, 0.0f); // Y axis - Green
        //axisAppearance.setColoringAttributes(ca);
        //ca.setColor(0.0f, 0.0f, 1.0f); // Z axis - Blue
    	
    	// X axis
        LineArray xAxis = new LineArray(2, LineArray.COORDINATES | LineArray.COLOR_3);
        xAxis.setCoordinate(0, new Point3f(0.0f, 0.0f, 0.0f));
        xAxis.setCoordinate(1, new Point3f(2.0f, 0.0f, 0.0f));
        xAxis.setColor(0, new Color3f(1.0f, 0.0f, 0.0f));
        xAxis.setColor(1, new Color3f(1.0f, 0.0f, 0.0f));
        Shape3D xAxisShape = new Shape3D(xAxis);
        objRoot.addChild(xAxisShape);

        // Y axis
        LineArray yAxis = new LineArray(2, LineArray.COORDINATES | LineArray.COLOR_3);
        yAxis.setCoordinate(0, new Point3f(0.0f, 0.0f, 0.0f));
        yAxis.setCoordinate(1, new Point3f(0.0f, 2.0f, 0.0f));
        yAxis.setColor(0, new Color3f(0.0f, 1.0f, 0.0f));
        yAxis.setColor(1, new Color3f(0.0f, 1.0f, 0.0f));
        Shape3D yAxisShape = new Shape3D(yAxis);
        objRoot.addChild(yAxisShape);

        // Z axis
        LineArray zAxis = new LineArray(2, LineArray.COORDINATES | LineArray.COLOR_3);
        zAxis.setCoordinate(0, new Point3f(0.0f, 0.0f, 0.0f));
        zAxis.setCoordinate(1, new Point3f(0.0f, 0.0f, 6.0f));
        zAxis.setColor(0, new Color3f(0.0f, 0.0f, 1.0f));
        zAxis.setColor(1, new Color3f(0.0f, 0.0f, 1.0f));
        Shape3D zAxisShape = new Shape3D(zAxis);
        objRoot.addChild(zAxisShape);
    }

    public static void main(String[] args) throws IOException {      
			new java3D_cylinder();
    }

}
