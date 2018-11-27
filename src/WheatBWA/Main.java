package WheatBWA;

import java.io.*;
import java.security.MessageDigest;
import java.security.DigestInputStream;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;


public class Main {

    public static void main(String[] args) {
        String filePath = "C:\\Users\\Van\\Desktop\\杂项";
        String outputPath = "C:\\Users\\Van\\Desktop\\md5.txt";

        try {
            File file = new File(outputPath);
            FileOutputStream out = new FileOutputStream(file);

            for (String i : FileListGenerator.getFileNameList(filePath, true)) {
                out.write((MD5Generator.fileMD5(filePath + "\\" + i) + " " + i + "\r")
                        .getBytes());
                System.out.println(MD5Generator.fileMD5(filePath + "\\" + i) + "  " + i);
            }
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}


/**
 * This is a class to generate MD5 value of string or file.
 * @author Duchangjiang Fan.
 **/
class MD5Generator {


    /**
     * This method is used to generate MD5 value of a string.
     * @param string : The string to get MD5 value.
     * @return : MD5 value in form of String.
     * */
    public static String stringMD5(String string) {

        try {
            MessageDigest md5_generator = MessageDigest.getInstance("MD5");

            md5_generator.update(string.getBytes());
            byte[] resultByteArray = md5_generator.digest();

            return MD5Generator.byteArrayToHex(resultByteArray);

        } catch (NoSuchAlgorithmException e) {
            e.printStackTrace();
            return "error";
        }

    }


    /**
     * This method is used to generate MD5 value of a file.
     * @param inputFile : Path of the file to get MD5 value.
     * @return : MD5 value in form of String.
     * */
    public static String fileMD5(String inputFile) {

        // A default read value is 256KB at once.
        int bufferSize = 256 * 1024;
        File file = new File(inputFile);
        FileInputStream fileInputStream = null;
        DigestInputStream digestInputStream = null;

        try {
            MessageDigest md5_generator = MessageDigest.getInstance("MD5");

            fileInputStream = new FileInputStream(file);
            digestInputStream = new DigestInputStream(fileInputStream, md5_generator);

            // Read the file.
            byte[] buffer = new byte[bufferSize];
            while (digestInputStream.read(buffer) > 0) ;

            // A similar process as MessageDigest.update()
            md5_generator = digestInputStream.getMessageDigest();

            // Store full digest into an byte array to transfer.
            byte[] resultByteArray = md5_generator.digest();

            return byteArrayToHex(resultByteArray);
        }

        catch (IOException | NoSuchAlgorithmException e) {
            e.printStackTrace();
            return null;
        }

        finally {

            try {
                digestInputStream.close();
                fileInputStream.close();
            }

            catch (NullPointerException | IOException e) {
                e.printStackTrace();
            }
        }
    }


    /**
     * This method is used to transform byte array to hexadecimal array, and a byte can be
     * transformed into two hexadecimal numbers.
     * @param byteArray : Byte code array to be transformed.
     * @return : Hexadecimal numbers in form of string.
     * */
    private static String byteArrayToHex(byte[] byteArray) {
        char[] hexDigits = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd',
                'e', 'f'};
        char[] resultCharArray = new char[byteArray.length * 2];

        int index = 0;
        for (byte b : byteArray) {
            resultCharArray[index++] = hexDigits[b >>> 4 & 0xf];
            resultCharArray[index++] = hexDigits[b & 0xf];
        }
        return new String(resultCharArray);
    }

}


/**
 * This is a class to obtain file list of a directory.
 * @author Duchangjiang Fan.
 **/
class FileListGenerator
{
    private static int stripLength;


    /**
     * This method is used to obtain file list of a directory recursively, so that files in
     * its subdirectory will be putted into the list too.
     * @param path : Path of the directory in form of string.
     * @param read_recursively : True for get files in subdirectory too.
     * @return : A list storing names of files(in form of xxx.xx), and files in subdirectory
     * recorded in a form of (xxx/xxx.xx).
     * */
    public static List<String> getFileNameList(String path, boolean read_recursively) {
        if (read_recursively) {
            if (stripLength == 0) {
                stripLength = path.length() + 1;
            }

            File file = new File(path);

            if (!file.exists()) {
                System.out.println(path + " not exists");
                return null;
            }
            if (file.isFile()) {
                System.out.println(path + " is a file");
                return null;
            }

            List<String> fileList = new ArrayList<>();
            File[] files = file.listFiles();

            for (File fs : files) {
                if (fs.isFile()) {
                    fileList.add(fs.getAbsolutePath().substring(stripLength));
                }
                else{
                    fileList.addAll(getFileNameList(fs.getPath(), true));
                }
            }
            return fileList;
        }

        return getFileNameList(path);
    }


    /**
     * This method is used to obtain file list of a directory excepting files in subdirectory.
     * @param path : Path of the directory in form of string.
     * @return : A list storing names of files(in form of xxx.xx).
     * */
    public static List<String> getFileNameList(String path) {

        File file = new File(path);

        if (!file.exists()) {
            System.out.println(path + " not exists");
            return null;
        }
        if (file.isFile()) {
            System.out.println(path + " is a file");
            return null;
        }

        File[] files = file.listFiles();

        List<String> fileList = new ArrayList<>();

        for (File fs : files) {
            if (fs.isFile()) {
                fileList.add(fs.getName());
            }
            else {
                fileList.add(fs.getName()+"\\");
            }
        }

        return fileList;
    }

}