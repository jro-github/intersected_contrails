package ContrailsWakeVortex2020;

import java.util.Scanner;

public class central_control_of_scenarios {
    /**
     * @param args
     */
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        int choice = -1;
        while (true) {
            System.out.println("Select scenario number to run:");
            System.out.println("===========================================================================");    
            System.out.println("0: Contrail 0 is independed with Contrail 1");
            System.out.println("1: Contrail 1 is created within Contrail 0");
            System.out.println("2: Contrail 0 and Contrail 1 is Parallel");
            System.out.println("3: Contrail 0 and Contrail 1 is Intersected with a set angle");
            System.out.println("===========================================================================");
            System.out.print("Please enter a scenario number (0-3): ");

            try {
                choice = Integer.parseInt(scanner.nextLine());
            } catch (NumberFormatException e) {
                System.out.println("===========================================================================");   
                System.out.println("Invalid input! Please enter a number between 0 and 3.");
                System.out.println("===========================================================================");   
                continue;
            }

            switch (choice) {
                case 0:
                    wake_vortex_individual_s0.run();
                    return;
                case 1:
                    wake_vortex_intersection_v2_for_s1.run();
                    return;
                case 2:
                    wake_vortex_intersection_v1_for_s2.run();
                    return;
                case 3:
                    double angle = 90; // default
                    while (true) {
                        System.out.print("Please enter the intersection angle in degrees (0-90, Default 90): ");
                        String input = scanner.nextLine();
                        if (input.trim().isEmpty()) {
                            // 90 Degrees by default
                            break;
                        }
                        try {
                            angle = Double.parseDouble(input);
                            if (angle >= 0 && angle <= 90) break;
                            System.out.println("Number must be between 0 and 90.");
                        } catch (NumberFormatException e) {
                            System.out.println("Invalid input! Please enter a number between 0 and 90.");
                        }
                    }
                    wake_vortex_intersection_v1_for_s3.run((int) angle);
                    return;
                default:
                    System.out.println("No such scenario!");
            }
        }
    }
}
