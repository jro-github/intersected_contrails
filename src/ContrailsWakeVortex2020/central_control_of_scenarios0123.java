package ContrailsWakeVortex2020;

public class central_control_of_scenarios0123 {
	public static void main(String[] args) {
		int choice = 1; // select to run which scenario.

        switch (choice) {
            case 0:
            	wake_vortex_individual_s0.run();
                break;
            case 1:
            	wake_vortex_intersection_v2_for_s1.run();
                break;
            case 2:
            	wake_vortex_intersection_v1_for_s2.run();
                break;
            case 3:
            	wake_vortex_intersection_v1_for_s3.run();
                break;
            default:
                System.out.println("No this scenario!");
        }
	}
}
