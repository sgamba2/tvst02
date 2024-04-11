
'''
from local_classes import *
# from mixing_inputs import *

class Project:
#------------------------------------------------------------------------------
# no need to have config files, can do initialization in python directly
#------------------------------------------------------------------------------
    def new_stage(self,name):
        self.fStage[name]            = Stage(name,self);
        return self.fStage[name]

    def dataset(self,dsid):
        return self.fDataset[dsid];
#------------------------------------------------------------------------------
# returns the name of the FCL file corresponding to the job - to be used by gen_fcl
#------------------------------------------------------------------------------
    def base_fcl(self,job,fcl_name):
        fmid = self.fFamilyID;              # familyID
        return self.fProjectName+'/datasets/'+fmid+'/'+job.stage().name()+'_'+fcl_name+'_'+fmid+'.fcl'

    def job_description(self,job):
        return self.fProjectName+'.'+job.input_dataset().id()+'.'+job.stage().name()+'_'+job.name()

    def __init__(self):

        project                      = 'tvst02'
        self.fFamilyID               = 'cals1b0'          # in fact, this is a family name
        self.fProjectName            = project;
        self.fStage                  = {}
        self.fDataset                = {};
        #------------------------------------------------------------------------------
        # datasets of this family
        # 1. stage 1 : generator input
        #------------------------------------------------------------------------------
        self.fDataset['cals1b0s00r0000'] = Dataset('generator'                   ,'cals1b0s00r0000','local')
        #------------------------------------------------------------------------------
        # 2. input for stage2 : datasets produced by stage1
        #------------------------------------------------------------------------------
        self.fDataset['cals1b0s11r0000'] = Dataset('sim.mu2e.cals1b0s11r0000.art','cals1b0s11r0000','local')
       
#------------------------------------------------------------------------------
# S1 10^10 cosmics
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s1');
        job                          = s.new_job('sim','cals1b0s00r0000');

        job.fRunNumber               = 1;
        job.fBaseFcl                 = self.base_fcl(job,'sim');

        job.fNInputFiles             = 1000                      # number of segments
                                     
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       = 10000000
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '40h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fMaxMemory               = '5000MB'

        odsid1                       = self.fFamilyID+'s11'+'r0000';
       

        job.fOutputStream            = ['PrimaryOutput'  ]
        job.fOutputDsID              = [odsid1 ] 
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0]]
        job.fOutputFormat            = ['art']
        
        # grid output dir
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;

# s2:dig : 
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s2');

        job                          = s.new_job('sim','cals1b0s00r0000');

        job.fBaseFcl                 = self.base_fcl(job,'sim');

        job.fNInputFiles             = -1                     # number of segments defined by s1:sim
             
        job.fMaxInputFilesPerSegment =  50
        job.fNEventsPerSegment       =  20000
        job.fResample                = 'no'   # yes/no        # for resampling, need to define the run number again
        job.fRequestedTime           = '3h'   
        job.fIfdh                    = 'xrootd'               # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        odsid21                      = self.fFamilyID+'s21'+'r0000';
        odsid22                      = self.fFamilyID+'s22'+'r0000';
        odsid23                      = self.fFamilyID+'s23'+'r0000';

        job.fOutputStream            = ['TargetStopOutput'            , 'ootStopOutput'               , 'IPAStopOutput'               ]
        job.fOutputDsID              = [odsid21                       , odsid22                       , odsid23                       ]
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0], 'sim.mu2e.'+job.fOutputDsID[1], 'sim.mu2e.'+job.fOutputDsID[2]]
        job.fOutputFormat            = ['art'                         , 'art'                         , 'art'                         ]

        # job description defined the grid output directory
        desc                         = project+'.'+job.input_dataset().id()+'.'+s.name()+'_'+job.name()
        job.fDescription             = desc;
 
        job                          = s.new_job('new_job','cals1b0s00r0000');

        job.fRunNumber               = 1;

        job.fNInputFiles             = 250#1000                      # number of segments
                                     
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       = 100000#10000000
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '10h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        odsid1                       = self.fFamilyID+'s12'+'r0000';
      

        job.fOutputStream            = ['PrimaryOutput'     ]
        job.fOutputDsID              = [odsid1                       ] 
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0]]
        job.fOutputFormat            = ['art'                          ]
#------------------------------------------------------------------------------
# init stage 2. a Stage can have one or several jobs associated with it
#------------------------------------------------------------------------------        
     s                            = self.new_stage('s2');
        job                          = s.new_job('sim','cals1b0s11r0000'); #'cals1b0s11r0000'

        job.fNInputFiles             = -1                     # number of segments defined by s1:sim
             
        job.fMaxInputFilesPerSegment =  1  #un evento per macchina
        job.fNEventsPerSegment       =  100000
        job.fResample                = 'no'   # yes/no        # for resampling, need to define the run number again
        job.fRequestedTime           = '10h'   
        job.fIfdh                    = 'xrootd'               # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        odsid21                      = self.fFamilyID+'s21'+'r0000';


        job.fOutputStream            = ['SignalOutput' ]
        job.fOutputDsID              = [odsid21        ]
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0]]
        job.fOutputFormat            = ['art'   ]

#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------
'''

#!/usr/bin/python

from local_classes import *
# from mixing_inputs import *

class Project(ProjectBase):

    def init_datasets(self):
#------------------------------------------------------------------------------
# datasets of this family
# 1. stage 1 : generator input
#------------------------------------------------------------------------------
        self.add_dataset(Dataset('generator'                          ,'cals1b0s00r0000','local'))
#------------------------------------------------------------------------------
# 2. input for stage2 : digit
#------------------------------------------------------------------------------
        self.add_dataset(Dataset('sim.mu2e.cals1b0s11r0000.tvst02.art','cals1b0s11r0000','local'))
        return

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
    def __init__(self,idsid=None):
        # print('Project.__init__: idsid=',idsid)
        ProjectBase.__init__(self,project='tvst02',family_id='cals1b0',idsid=idsid);
        self.init_datasets()
#------------------------------------------------------------------------------
# S1 10^8 proton interactions in the PT, half field in the DS
#------------------------------------------------------------------------------        
        s                            = self.new_stage('s1');
        job                          = s.new_job('sim','cals1b0s00r0000');

        job.fRunNumber               = 1;

        job.fNInputFiles             = 1000#1000                      # number of segments
                                     
        job.fMaxInputFilesPerSegment =  1
        job.fNEventsPerSegment       = 2000000#10000000
        job.fResample                = 'no'   # yes/no
        job.fRequestedTime           = '40h'
        job.fIfdh                    = 'xrootd'                 # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        odsid1                       = self.fFamilyID+'s11'+'r0000';
      

        job.fOutputStream            = ['PrimaryOutput'     ]
        job.fOutputDsID              = [odsid1                       ] 
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0]]
        job.fOutputFormat            = ['art'                          ]
#



 
#------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------
        s                            = self.new_stage('s2');
        job                          = s.new_job('sim','cals1b0s11r0000'); #'cals1b0s11r0000'

        job.fNInputFiles             = -1                     # number of segments defined by s1:sim
             
        job.fMaxInputFilesPerSegment =  1  #un evento per macchina
        job.fNEventsPerSegment       =  100000
        job.fResample                = 'no'   # yes/no        # for resampling, need to define the run number again
        job.fRequestedTime           = '30h'   
        job.fIfdh                    = 'xrootd'               # ifdh/xrootd
        job.fMaxMemory               = '3000MB'

        odsid21                      = self.fFamilyID+'s21'+'r0000';


        job.fOutputStream            = ['TCalStat' ]
        job.fOutputDsID              = [odsid21        ]
        job.fOutputFnPattern         = ['sim.mu2e.'+job.fOutputDsID[0]]
        job.fOutputFormat            = ['art'   ]
