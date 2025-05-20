#include "strength.h"
#include "unity.h"
#include <stdlib.h>

void test_calculate_equivalent_strain(void) {
    // Test case values
    // These values are from a real-world example of a uniaxial tensile test
    double strain[6] = {0.1, -0.05, -0.05, 0.0, 0.0, 0.0};
    
    // Expected result (calculated manually or from a verified source)
    double expected_eq_strain = 0.1732050807568877;  // √(2/3 * (0.1 - (-0.05))^2 + (-0.05 - (-0.05))^2 + (-0.05 - 0.1)^2)
    
    // Calculate the result using our function
    double calculated_eq_strain = equiv_strain(strain);
    
    // Check if the calculated result matches the expected result
    TEST_ASSERT_FLOAT_WITHIN(1e-10, expected_eq_strain, calculated_eq_strain);
}

void test_calculate_equivalent_strain_rate(void) {
    // Test case values
    // These values represent a simple shear deformation scenario
    double strain_rate[6] = {0.0, 0.0, 0.0, 0.1, 0.0, 0.0};
    
    // Expected result (calculated manually or from a verified source)
    // For simple shear, the equivalent strain rate is γ̇ / √3
    double expected_eq_strain_rate = 0.1 / sqrt(3);
    
    // Calculate the result using our function
    double calculated_eq_strain_rate = equiv_strain(strain_rate);
    
    // Check if the calculated result matches the expected result
    TEST_ASSERT_FLOAT_WITHIN(1e-10, expected_eq_strain_rate, calculated_eq_strain_rate);
}

void test_calculate_equivalent_strain_rate_uniaxial(void) {
    // Test case for uniaxial tension/compression
    double strain_rate[6] = {0.2, -0.1, -0.1, 0.0, 0.0, 0.0};
    
    // Expected result
    // For uniaxial deformation, the equivalent strain rate is equal to the absolute value of the axial strain rate
    double expected_eq_strain_rate = 0.2;
    
    // Calculate the result using our function
    double calculated_eq_strain_rate = equiv_strain(strain_rate);
    
    // Check if the calculated result matches the expected result
    TEST_ASSERT_FLOAT_WITHIN(1e-10, expected_eq_strain_rate, calculated_eq_strain_rate);
}

int main() {
    UNITY_BEGIN();
    RUN_TEST(test_calculate_equivalent_strain);
    RUN_TEST(test_calculate_equivalent_strain_rate);
    RUN_TEST(test_calculate_equivalent_strain_rate_uniaxial);
    return UNITY_END();
}