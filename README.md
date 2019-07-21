# PrivatePlantGenetics
Some easy code  written by myself ~

# My Release 1: 
## Aoyue-monitor 
[https://github.com/AoyueBi/PrivateJavaCodeByAoyue/releases](https://github.com/AoyueBi/PrivateJavaCodeByAoyue/releases)

**Hello!** To use this **Aoyue-monitor** jar, please add **4** parameters and prepare **two** files.

* @param1: keepjobnum eg: **20** 
* @param2: sleeptime eg: **5** 
* @param3: getJobRunshS eg: **/Users/Aoyue/Documents/getJobRun.sh**
* @param4: allshS eg: **/Users/Aoyue/Documents/001\_allsh_samtoolsView.sh**

### Usage
From linux command line: 
`java -jar Aoyue-monitor.jar 20 5 /pathto/getJobRun.sh /pathto/001_allsh_samtoolsView.sh > log_monitor.txt &`

### Note

- **@param1**: The job you wanna keep in linux process 

- **@param2**: Unit: min, the interval time you wanna check the ps aux process
- **@param3**: The sh file to produce the job num running now eg:`ps aux|grep aoyue|grep samtools|grep view|wc -l `
- **@param4**: All the CMDs in one script, one line one sample, please DO NOT add & symbol at the end of each line.

### Work Flowchart

![Flowchart]
(https://ibb.co/vB67G3m)

### Contributor
Aoyue Bi

[biaoyue17@genetics.ac.cn]()

[https://plantgeneticslab.weebly.com/](https://plantgeneticslab.weebly.com/)

[https://aoyuebi.wixsite.com/bay9624](https://aoyuebi.wixsite.com/bay9624)

Thanks Dr. Gao (Qiang Gao in BGI) for helping me organizing my ideas.
