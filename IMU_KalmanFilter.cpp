/**
 *Bibliothek für eine Sammlung von IMU Funktionen nur Bestimmung der Lage aus Sensordaten
 *
 * Stefan Schlee, 2020
 */

#include "IMU_KalmanFilter.h"

/**
 * Erstellt eine IMU Instanz und parametriert den Filter
 * Q..: Varianzen der Zustände
 * R..: Varianzen der Messungen
 */
Attitude_3D_Kalman::Attitude_3D_Kalman(float Abtastzeit_s, float Qyaw, float Qpitch_roll, float Qgyrobias, float Ryaw, float Rpitch_roll)
    : Abtastzeit_s(Abtastzeit_s), Qyaw(Qyaw*Abtastzeit_s), Qpitch_roll(Qpitch_roll*Abtastzeit_s), Qgyrobias(Qgyrobias*Abtastzeit_s), Ryaw(Ryaw), Rpitch_roll(Rpitch_roll)
{
}

/**
 * Kalman Filter Algorithmus
 * Acc: Beschleunigung in g
 * Gyro: Winkelgeschwindigkeit in degrees per second
 * Mag: Magnetometer, Skalierung egal
 * 
 * Koordinatensystem wie das ungedrehte Roll pitch yaw System (x y z)
 */
Attitude_3D_t Attitude_3D_Kalman::update(float Acc_X, float Acc_Y, float Acc_Z, float Gyro_X, float Gyro_Y, float Gyro_Z, float Mag_X, float Mag_Y, float Mag_Z)
{

    /*********  Variablen  **********/
    //Neue Messung, yaw pitch roll
    float measurement[3];

    //Innovation, yaw pitch roll
    float y[3];

    //Kalman Gain, Zustand x Innovation
    float K1_1, K6_1, K2_2, K5_2, K6_2, K3_3, K4_3, K5_3, K6_3;

    /*********  Berechnung der neuen Messungen  **********/
    //Winkel in RAD für Winkelfunktionen
    float pitch_rad = pitch * DEG_TO_RAD;
    float roll_rad = roll * DEG_TO_RAD;

    //Pitch und Roll aus Accelerometer
    measurement[2] = atanf(Acc_Y / sqrtf(Acc_X * Acc_X + Acc_Z * Acc_Z)) * RAD_TO_DEG;
    measurement[1] = atanf(-1 * Acc_X / sqrtf(Acc_Y * Acc_Y + Acc_Z * Acc_Z)) * RAD_TO_DEG;

    //Magnetometer Messung in horizontale Ebene transformieren
    float XH = Mag_X * cosf(pitch_rad) + Mag_Z * cosf(roll_rad) * sinf(pitch_rad) + Mag_Y * sinf(pitch_rad) * sinf(roll_rad);
    float YH = Mag_Y * cosf(roll_rad) - Mag_Z * sinf(roll_rad);

    //Yaw mit Magnetometer berechnen
    measurement[0] = atan2f(XH, YH) * RAD_TO_DEG;

    //prevent 360deg jump
    if ((measurement[0] - yaw) > 180)
    {
        measurement[0] -= 360;
    }
    else if ((measurement[0] - yaw) < -180)
    {
        measurement[0] += 360;
    }

    /*********  Kalman Algorithmus  **********/
    //new state estimates mit Vorwärtstransitions Matrix
    yaw += ((cosf(roll_rad) * (Gyro_Z - zbias)) / cosf(pitch_rad) + (sinf(roll_rad) * (Gyro_Y - ybias)) / cosf(pitch_rad)) * Abtastzeit_s;
    pitch += (cosf(roll_rad) * (Gyro_Y - ybias) - sinf(roll_rad) * (Gyro_Z - zbias)) * Abtastzeit_s;
    roll += (Gyro_X - xbias + cosf(roll_rad) * tanf(pitch_rad) * (Gyro_Z - zbias) + tanf(pitch_rad) * sinf(roll_rad) * (Gyro_Y - ybias)) * Abtastzeit_s;

    //new error covariance estimate
    float P4_4_temp = P4_4;
    float P5_5_temp = P5_5;
    float P6_6_temp = P6_6;

    P1_1 = P1_1 + Qyaw - (P1_6 * Abtastzeit_s * cosf(roll_rad)) / cosf(pitch_rad) - (Abtastzeit_s * cosf(roll_rad) * (P1_6 - (P6_6 * Abtastzeit_s * cosf(roll_rad)) / cosf(pitch_rad))) / cosf(pitch_rad);

    P2_2 = P2_2 + Qpitch_roll - Abtastzeit_s * cosf(roll_rad) * (P2_5 - P5_5 * Abtastzeit_s * cosf(roll_rad)) + Abtastzeit_s * sinf(roll_rad) * (P2_6 + P6_6 * Abtastzeit_s * sinf(roll_rad)) - P2_5 * Abtastzeit_s * cosf(roll_rad) + P2_6 * Abtastzeit_s * sinf(roll_rad);

    P3_3 = P3_3 + Qpitch_roll - P3_4 * Abtastzeit_s - Abtastzeit_s * (P3_4 - P4_4 * Abtastzeit_s) - Abtastzeit_s * tanf(pitch_rad) * sinf(roll_rad) * (P3_5 - P5_5 * Abtastzeit_s * tanf(pitch_rad) * sinf(roll_rad)) - P3_6 * Abtastzeit_s * cosf(roll_rad) * tanf(pitch_rad) - P3_5 * Abtastzeit_s * tanf(pitch_rad) * sinf(roll_rad) - Abtastzeit_s * cosf(roll_rad) * tanf(pitch_rad) * (P3_6 - P6_6 * Abtastzeit_s * cosf(roll_rad) * tanf(pitch_rad));

    P4_4 += Qgyrobias;
    P5_5 += Qgyrobias;
    P6_6 += Qgyrobias;

    P1_6 -= (P6_6_temp * Abtastzeit_s * cosf(roll_rad)) / cosf(pitch_rad);
    P2_5 -= P5_5_temp * Abtastzeit_s * cosf(roll_rad);
    P2_6 += P6_6_temp * Abtastzeit_s * sinf(roll_rad);
    P3_4 -= P4_4_temp * Abtastzeit_s;
    P3_5 -= P5_5_temp * Abtastzeit_s * tanf(pitch_rad) * sinf(roll_rad);
    P3_6 -= P6_6_temp * Abtastzeit_s * cosf(roll_rad) * tanf(pitch_rad);

    //Innovation berechnen
    y[0] = measurement[0] - yaw;
    y[1] = measurement[1] - pitch;
    y[2] = measurement[2] - roll;

    //Kalman Gain berechnen
    K1_1 = P1_1 / (P1_1 + Ryaw);
    K2_2 = P2_2 / (P2_2 + Rpitch_roll);
    K3_3 = P3_3 / (P3_3 + Rpitch_roll);
    K4_3 = P3_4 / (P3_3 + Rpitch_roll);
    K5_2 = P2_5 / (P2_2 + Rpitch_roll);
    K5_3 = P3_5 / (P3_3 + Rpitch_roll);
    K6_1 = P1_6 / (P1_1 + Ryaw);
    K6_2 = P2_6 / (P2_2 + Rpitch_roll);
    K6_3 = P3_6 / (P3_3 + Rpitch_roll);

    //corrected state estimate
    yaw += K1_1 * y[0];
    pitch += K2_2 * y[1];
    roll += K3_3 * y[2];

    xbias += K4_3 * y[2];
    ybias += K5_2 * y[1] + K5_3 * y[2];
    zbias += K6_1 * y[0] + K6_2 * y[1] + K6_3 * y[2];

    //corrected error covariance matrix
    P1_1 -= P1_1 * K1_1;
    P2_2 -= P2_2 * K2_2;
    P3_3 -= P3_3 * K3_3;
    P4_4 -= P3_4 * K4_3;
    P5_5 -= (K5_2 * P2_5 + K5_3 * P3_5);
    P6_6 -= (K6_1 * P1_6 + K6_2 * P2_6 + K6_3 * P3_6);

    P1_6 -= P1_6 * K1_1;
    P2_5 -= P2_5 * K2_2;
    P2_6 -= P2_6 * K2_2;
    P3_4 -= P3_4 * K3_3;
    P3_5 -= P3_5 * K3_3;
    P3_6 -= P3_6 * K3_3;

    //contrain yaw to +-180deg
    if (yaw > 180)
    {
        yaw -= 360;
    }
    else if (yaw < -180)
    {
        yaw += 360;
    }

    //Output Struktur erzeugen und zurückgeben
    return (Attitude_3D_t){yaw, pitch, roll, Gyro_X-xbias, Gyro_Y-ybias, Gyro_Z-zbias};
}