/**
 *Bibliothek für eine Sammlung von IMU Funktionen nur Bestimmung der Lage aus Sensordaten
 *
 * Stefan Schlee, 2020
 */
#ifndef IMU_KALMANFILTER_H
#define IMU_KALMANFILTER_H
#include <Arduino.h>

/**
 * Return Daten eine 3D Attitude IMUs
 * roll nach vorne, yaw nach oben
 */
typedef struct Attitude_3D_t
{
    float yaw;
    float pitch;
    float roll;
    float unbiased_gyro_x;
    float unbiased_gyro_y;
    float unbiased_gyro_z;
} Attitude_3D_t;

/**
 * 3D_Attitude IMU mit einem Kalman Filter
 */
class Attitude_3D_Kalman
{
public:
    Attitude_3D_Kalman(float Abtastzeit_s, float Qyaw, float Qpitch_roll, float Qgyrobias, float Ryaw, float Rpitch_roll);
    Attitude_3D_t update(float Acc_X, float Acc_Y, float Acc_Z, float Gyro_X, float Gyro_Y, float Gyro_Z, float Mag_X, float Mag_Y, float Mag_Z);

private:
    //Abtastzeit des Filters für Integration [Sekunden]#
    const float Abtastzeit_s;

    //Varianzen der Zustände
    const float Qyaw, Qpitch_roll, Qgyrobias; 

    //Varianzen der Messungen
    const float Ryaw, Rpitch_roll;            

    //Zustände
    float yaw;
    float pitch;
    float roll;
    float xbias;
    float ybias;
    float zbias;

    //error covariance matrix, nur elemente ungleich 0, außerdem symmetrische Matrix
    float P1_1, P2_2, P3_3, P4_4, P5_5, P6_6;
    float P3_4, P3_5, P3_6, P2_5, P2_6, P1_6;
};

#endif // IMU_KALMANFILTER_H
