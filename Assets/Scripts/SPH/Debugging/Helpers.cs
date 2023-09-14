using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using SerializableTypes;
using Random = System.Random;
using System.Linq;

namespace Helpers {
    [System.Serializable]
    public static class HelperMethods {
        //public static Random rng = new Random();  

        public static bool HasComponent<T> (GameObject obj, out T toReturn) {
            toReturn = obj.GetComponent<T>();
            return toReturn != null;
        }
        public static bool HasComponent<T> (GameObject obj) {
            return obj.GetComponent<T>() != null;
        }
        public static bool HasComponent<T> (Transform t, out T toReturn) {
            return HasComponent<T>(t.gameObject, out toReturn);
        }
        public static bool HasComponent<T> (Transform t) {
            return HasComponent<T>(t.gameObject);
        }

        public static float Map(float value, float from1, float to1, float from2, float to2) {
            return (value - from1) / (to1 - from1) * (to2 - from2) + from2;
        }

        // Adapted from Fisher-Yates Shuffle, code adopted from https://stackoverflow.com/questions/273313/randomize-a-listt
        public static List<T> Shuffle<T>(this List<T> list) {
            Random rng = new Random();
            List<T> newList = new List<T>(list);
            int n = newList.Count;
            while (n > 1) {
                n--; 
                int k = rng.Next(n + 1); 
                T value = newList[k]; 
                newList[k] = newList[n]; 
                newList[n] = value; 
            }
            return newList;
        }
        public static List<T> Flatten2D<T>(this List<List<T>> list) {
            List<T> flatten = new List<T>();
            foreach(List<T> nested in list) {
                foreach(T item in nested) {
                    flatten.Add(item);
                }
            }
            return flatten;
        }
        public static List<T2> Flatten2D<T1,T2>(this Dictionary<T1, List<T2>> dict) {
            List<T2> flatten = new List<T2>();
            foreach(KeyValuePair<T1, List<T2>> kvp in dict) {
                foreach(T2 item in kvp.Value) {
                    flatten.Add(item);
                }
            }
            return flatten;
        }

        // Adapted from: https://stackoverflow.com/questions/19141259/how-to-enqueue-a-list-of-items-in-c
        public static void AddRange<T>(this Queue<T> queue, IEnumerable<T> enu) {
            foreach (T obj in enu)
                queue.Enqueue(obj);
        }

        public static bool Compare(this SVector3 original, SVector3 other) {
            return original.x == other.x && original.y == other.y && original.z == other.z;
        }
        public static bool Compare(this SVector3 original, Vector3 other) {
            return original.x == other.x && original.y == other.y && original.z == other.z;
        }
        public static bool Compare(this SVector4 original, SVector4 other) {
            return original.x == other.x && original.y == other.y && original.z == other.z && original.w == other.w;
        }
        public static bool Compare(this SVector4 original, Vector4 other) {
            return original.x == other.x && original.y == other.y && original.z == other.z && original.w == other.w;
        }

        // Adapted from: https://answers.unity.com/questions/288338/how-do-i-compare-quaternions.html
        public static bool Compare(this Quaternion original, Quaternion other, float acceptableRange) {
            return 1 - Mathf.Abs(Quaternion.Dot(original, other)) < acceptableRange;
        }

        // Adapted from: https://stackoverflow.com/questions/1792470/subset-of-array-in-c-sharp
        public static T[] RangeSubset<T>(this T[] array, int startIndex, int length) {
            T[] subset = new T[length];
            Array.Copy(array, startIndex, subset, 0, length);
            return subset;
        }

        // [Source: https://forum.unity.com/threads/encoding-vector2-and-vector3-variables-into-single-int-or-float-and-back.448346/]
        // [Commentor: emotitron]
        public static int encodeVector3ToInt(Vector3 v) {
            //Vectors must stay within the -512 to 512 range per axis - no error handling coded here
            //Add 512 to get numbers into the 0-1024 range rather than -512 to 512 range
            //Multiply by 10 to save one decimal place from rounding
            int xcomp = Mathf.RoundToInt((v.x * 10)) + 512;
            int ycomp = Mathf.RoundToInt((v.y * 10)) + 512;
            int zcomp = Mathf.RoundToInt((v.z * 10)) + 512;
            return xcomp + ycomp * 1024 + zcomp * 1048576;
        }
        public static Vector3 decodeVector3FromInt(int i) {
            //Get the leftmost bits first. The fractional remains are the bits to the right.
            // 1024 is 2 ^ 10 - 1048576 is 2 ^ 20 - just saving some calculation time doing that in advance
            float z = Mathf.Floor(i / 1048576);
            float y = Mathf.Floor ((i - z * 1048576) / 1024);
            float x = (i - y * 1024 - z * 1048576);
            // subtract 512 to move numbers back into the -512 to 512 range rather than 0 - 1024
            return new Vector3 ((x - 512) / 10, (y - 512) / 10, (z - 512) / 10);
        }
        public static SVector3 decodeSVector3FromInt(int i) {
            //Get the leftmost bits first. The fractional remains are the bits to the right.
            // 1024 is 2 ^ 10 - 1048576 is 2 ^ 20 - just saving some calculation time doing that in advance
            float z = Mathf.Floor(i / 1048576);
            float y = Mathf.Floor ((i - z * 1048576) / 1024);
            float x = (i - y * 1024 - z * 1048576);
            // subtract 512 to move numbers back into the -512 to 512 range rather than 0 - 1024
            return new SVector3 ((x - 512) / 10, (y - 512) / 10, (z - 512) / 10);
        }
    }

    [System.Serializable]
    public class SaveSystemMethods {
        public static string GetSaveLoadDirectory(string path = "") {
            return (path != null && path.Length > 0) 
                ? (path.EndsWith("/")) 
                    ? Application.dataPath + "/" + path 
                    : Application.dataPath + "/" + path + "/" 
                : Application.dataPath + "/";
        }
        public static bool CheckDirectoryExists(string dirPath) {
            return Directory.Exists(dirPath);
        }
        public static bool CheckOrCreateDirectory(string dirPath) {
            if (!CheckDirectoryExists(dirPath)) Directory.CreateDirectory(dirPath);
            return true;
        }
        public static bool CheckFileExists(string filePath) {
            return File.Exists(filePath);
        }
        public static bool DeleteDirectory(string dirPath) {
            if (CheckDirectoryExists(dirPath)) {
                Directory.Delete(dirPath,true);
                return true;
            }
            return false;
        }
        public static bool DeleteFile(string filePath) {
            if (CheckFileExists(filePath)) {
                File.Delete(filePath);
                return true;
            }
            return false;
        }
        
        public static string ConvertToJSON<T>(T data) {
            return JsonUtility.ToJson(data, true);
        }
        public static T ConvertFromJSON<T>(string data) {
            return JsonUtility.FromJson<T>(data);
        }
        public static bool SaveJSON(string filePath, string json) {
            if (filePath.EndsWith(".json")) File.WriteAllText(filePath, json);
            else File.WriteAllText(filePath + ".json", json);
            return true;
        }
        public static bool LoadJSON<T>(string filePath, out T output) {
            string actualFilePath = (filePath.EndsWith(".json")) ? filePath : filePath+".json";
            if (CheckFileExists(actualFilePath)) {
                string fileContents = File.ReadAllText(actualFilePath);
                output = ConvertFromJSON<T>(fileContents);
                return true;
            } else {
                output = default(T);
                return false;
            }
        }
        public static bool SaveCSV<T>(string filePath, List<string> header, List<T> data) {
            string p = (filePath.EndsWith(".csv")) ? filePath : filePath+".csv";
            StreamWriter writer = new StreamWriter(p);
            writer.WriteLine(string.Join(",", header));
            foreach(T line in data) {
                List<string> lineContent = new List<string>();
                foreach(string headerVal in header) {
                    lineContent.Add(line.GetType().GetField(headerVal).GetValue(line).ToString());
                }
                writer.WriteLine(string.Join(",",lineContent));
            }
            writer.Flush();
            writer.Close();
            return true;
        }
        
        public static string[] ReadCSV(TextAsset asset) {
            return asset.text.Split(new string[] {",","\n"}, StringSplitOptions.None);
        }
    }

    // Source: https://forum.unity.com/threads/debug-drawarrow.85980/
    [System.Serializable]
    public class DrawingMethods {
        public static void ForGizmos(Vector3 pos, Vector3 direction, float arrowHeadLength = 0.25f, float arrowHeadAngle = 20.0f) {
            Gizmos.DrawRay(pos, direction);

            Vector3 right = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180+arrowHeadAngle,0) * new Vector3(0,0,1);
            Vector3 left = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180-arrowHeadAngle,0) * new Vector3(0,0,1);
            Gizmos.DrawRay(pos + direction, right * arrowHeadLength);
            Gizmos.DrawRay(pos + direction, left * arrowHeadLength);
        }
 
        public static void ForGizmos(Vector3 pos, Vector3 direction, Color color, float arrowHeadLength = 0.25f, float arrowHeadAngle = 20.0f) {
            Gizmos.color = color;
            Gizmos.DrawRay(pos, direction);
       
            Vector3 right = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180+arrowHeadAngle,0) * new Vector3(0,0,1);
            Vector3 left = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180-arrowHeadAngle,0) * new Vector3(0,0,1);
            Gizmos.DrawRay(pos + direction, right * arrowHeadLength);
            Gizmos.DrawRay(pos + direction, left * arrowHeadLength);
        }
 
        public static void ForDebug(Vector3 pos, Vector3 direction, float arrowHeadLength = 0.25f, float arrowHeadAngle = 20.0f) {
            Debug.DrawRay(pos, direction);
       
            Vector3 right = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180+arrowHeadAngle,0) * new Vector3(0,0,1);
            Vector3 left = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180-arrowHeadAngle,0) * new Vector3(0,0,1);
            Debug.DrawRay(pos + direction, right * arrowHeadLength);
            Debug.DrawRay(pos + direction, left * arrowHeadLength);
        }
    
        public static void ForDebug(Vector3 pos, Vector3 direction, Color color, float arrowHeadLength = 0.25f, float arrowHeadAngle = 20.0f) {
            Debug.DrawRay(pos, direction, color);
       
            Vector3 right = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180+arrowHeadAngle,0) * new Vector3(0,0,1);
            Vector3 left = Quaternion.LookRotation(direction) * Quaternion.Euler(0,180-arrowHeadAngle,0) * new Vector3(0,0,1);
            Debug.DrawRay(pos + direction, right * arrowHeadLength, color);
            Debug.DrawRay(pos + direction, left * arrowHeadLength, color);
        }
    }
}