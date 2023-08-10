using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Networking;

namespace WebTools {
    [System.Serializable]
    public static class WebRequests {
        
        public class Response {
            public bool success = false;
            public string response = null;
            public Response(bool success, string response=null) {
                this.success = success;
                this.response = response;
            } 
        }

        public static IEnumerator GetRequest(string url, bool verbose=false, Action<Response> callback = null) {
            Response r;
            UnityWebRequest uwr = UnityWebRequest.Get(url);
            yield return uwr.SendWebRequest();

            if (uwr.result == UnityWebRequest.Result.ConnectionError) {
                // Error has occurred. Prints error message and returns false
                r = new Response(false, uwr.error.ToString());
                if (verbose) Debug.LogError($"Error: GET Request to {url} returned error: {uwr.error}");
            } else {
                // If we made it through, we pass true
                r = new Response(true, uwr.downloadHandler.text);
                if (verbose) Debug.Log($"GET Request to {url} Success: {uwr.downloadHandler.text}");
            }
            
            uwr.Dispose();
            if (callback != null) callback(r);
        }

        public static IEnumerator PostRequestWithForm(string url, Dictionary<string,string> data, bool verbose=false, Action<Response> callback = null) {
            Response r;

            WWWForm form = new WWWForm();
            foreach(KeyValuePair<string,string> kvp in data) {
                form.AddField(kvp.Key, kvp.Value);
            }
            
            UnityWebRequest uwr = UnityWebRequest.Post(url, form);
            yield return uwr.SendWebRequest();

            if (uwr.result == UnityWebRequest.Result.ConnectionError) {
                // Error has occurred. Prints error message and returns false
                r = new Response(false, uwr.error.ToString());
                if (verbose) Debug.LogError($"Error: POST Request to {url} returned error: {uwr.error}");
            } else {
                // If we made it through, we pass true
                r = new Response(true, uwr.downloadHandler.text);
                if (verbose) Debug.Log($"GET Request to {url} Success: {uwr.downloadHandler.text}");
            }
            
            uwr.Dispose();
            if (callback != null) callback(r);
        }

        public static IEnumerator PostRequestWithJSON(string url, string data, bool verbose = false, Action<Response> callback = null) {
            Response r;
            var uwr = new UnityWebRequest(url, "POST");
            byte[] jsonToSend = new System.Text.UTF8Encoding().GetBytes(data);
            uwr.uploadHandler = (UploadHandler)new UploadHandlerRaw(jsonToSend);
            uwr.downloadHandler = (DownloadHandler)new DownloadHandlerBuffer();
            uwr.SetRequestHeader("Content-Type", "application/json");

            //Send the request then wait here until it returns
            yield return uwr.SendWebRequest();
            
            if (uwr.result == UnityWebRequest.Result.ConnectionError) {
                // Error has occurred. Prints error message and returns false
                r = new Response(false, uwr.error.ToString());
                if (verbose) Debug.LogError($"Error: POST Request to {url} returned error: {uwr.error}");
            } else {
                // If we made it through, we pass true
                r = new Response(true, uwr.downloadHandler.text);
                if (verbose) Debug.Log($"GET Request to {url} Success: {uwr.downloadHandler.text}");
            }

            uwr.Dispose();
            if (callback != null) callback(r);
        }
    }   
}