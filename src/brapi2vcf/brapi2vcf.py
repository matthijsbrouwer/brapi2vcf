"""Main BrAPI2vcf module"""
import logging
import sys
import requests
import tempfile
import shutil
from datetime import datetime
from contextlib import nullcontext

class BrAPI2vcf:
    """Main BrAPI2vcf class"""

    def __init__(self,url):
        self._url = str(url)
        #set logging
        self._logger = logging.getLogger("brapi2vcf")
        #check serverinfo
        self._serverinfo = self._brapiRequest("serverinfo")


    def _brapiRequest(self,call,params={}):
        headers = {"Accept": "application/json"}
        url = "%s/%s" % (self._url,call)
        response = requests.get(url, params=params, headers=headers)
        return response.json()
    
    def vcf(self, outputFile: str = None):
        buffer = tempfile.TemporaryFile(mode="w+")
        try:
            hasAllelematrix = False
            hasVariants = False
            for callEntry in self._serverinfo.get("result",{}).get("calls",[]):
                if callEntry["service"]=="allelematrix":
                    hasAllelematrix = True
                elif callEntry["service"]=="variants":
                    hasVariants = True
            if hasAllelematrix:
                with open(outputFile, "w") if not outputFile is None else nullcontext(sys.stdout) as output:
                    vcf_format_lines = {}
                    vcf_info_lines = {
                        "svtype": "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                        "end": "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
                        "svlen": "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
                        "cipos": "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">",
                        "ciend": "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">"
                    }
                    output.write("##fileformat=VCFv4.3\n")
                    now = datetime.now()
                    output.write("##fileDate=%s\n" % now.strftime("%Y%m%d"))
                    output.write("##source=%s\n" % self._serverinfo.get("result",{}).get("serverName","BrAPI"))
                    variants = {}
                    sampleDbIds = []
                    #get samples from allelematrix
                    page = 0
                    totalPages = 1
                    while page<totalPages:
                        response = self._brapiRequest("allelematrix",{"dimensionVariantPageSize": 1, 
                                                                      "dimensionCallSetPage": page,
                                                                      "dimensionCallSetPageSize": 1000, 
                                                                      "preview": True})
                        for entry in response.get("result",{}).get("pagination",[]):
                            if entry.get("dimension","")=="CALLSETS":
                                totalPages = entry.get("totalPages",1)
                        sampleDbIds.extend(response.get("result",{}).get("callSetDbIds",[]))
                        page+=1
                    buffer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
                    for sampleDbId in sampleDbIds:
                        buffer.write("\t%s" % str(sampleDbId))
                    buffer.write("\n")
                    #get variant info
                    if hasVariants:
                        page = 0
                        totalPages = 1
                        while page<totalPages:
                            response = self._brapiRequest("variants",{"page": page})
                            totalPages = response.get("metadata",{}).get("pagination",{}).get("totalPages",1)
                            for entry in response.get("result",{}).get("data",[]):
                                if "variantDbId" in entry:
                                    variantOutputLine = variantOutputLine = [
                                        ".",
                                        ".",
                                        entry["variantDbId"],
                                        ".",
                                        ".",
                                        ".",
                                        "",
                                        "",
                                        []]
                                    start = entry.get("start",None)
                                    if entry.get("referenceName",None):
                                        variantOutputLine[0] = str(entry.get("referenceName","")).strip()
                                    elif entry.get("referenceDbId",None):
                                        variantOutputLine[0] = str(entry.get("referenceDbId")).strip()
                                    if isinstance(entry.get("start",None),int):
                                        variantOutputLine[1] = str(entry["start"])
                                    if entry.get("referenceBases",None):
                                        variantOutputLine[3] = str(entry.get("referenceBases")).strip()
                                    if len(entry.get("alternateBases",[]))>0:
                                        variantOutputLine[4] = ",".join([str(item) for item in entry.get("alternateBases")])
                                    if entry.get("filtersApplied",False):
                                        if entry.get("filtersPassed",False):
                                            variantOutputLine[6] = "PASS"
                                        elif len(entry.get("filtersFailed",[]))>0:
                                            variantOutputLine[4] = ";".join([str(item) for item in entry.get("filtersFailed")])
                                    else:
                                        variantOutputLine[6] = "."
                                    infoList = []
                                    if isinstance(entry.get("end",None),int):
                                        infoList.append("END=%s" % str(entry["end"]))
                                    if isinstance(entry.get("svlen",None),int):
                                        infoList.append("SVLEN=%s" % str(entry["svlen"]))
                                    if len(entry.get("cipos",[]))>0:
                                        infoList.append("CIPOS=%s" %",".join([str(item) for item in entry.get("cipos")]))
                                    if len(entry.get("ciend",[]))>0:
                                        infoList.append("CIEND=%s" % ",".join([str(item) for item in entry.get("ciend")]))
                                    if entry.get("variantType",None):
                                        infoList.append("SVTYPE=%s" % str(entry["variantType"]))
                                    variantOutputLine[7] = ";".join(infoList)
                                    variants[entry["variantDbId"]] = variantOutputLine
                            page+=1
                    #get allelematrix
                    variantPage = 0
                    totalVariantPages = 1
                    totalCallSetPages = 1
                    while variantPage<totalVariantPages:
                        callSetPage = 0
                        variantOutputLines = []
                        sampleCounterOffset = 9
                        while callSetPage<totalCallSetPages:
                            params = {"dimensionVariantPage": variantPage, "dimensionVariantPageSize":100, 
                                      "dimensionCallSetPage": callSetPage, "dimensionCallSetPageSize":100,
                                      "sepPhased": "/", "sepUnphased": "|", "unknownString": "."}
                            response = self._brapiRequest("allelematrix",params)
                            #update pagination
                            for entry in response.get("result",{}).get("pagination",[]):
                                if entry.get("dimension","")=="VARIANTS":
                                    totalVariantPages = entry.get("totalPages",1)
                                elif entry.get("dimension","")=="CALLSETS":
                                    totalCallSetPages = entry.get("totalPages",1)
                            #output variant id
                            if callSetPage==0:
                                variantOutputLines = []                                
                                for variantDbId in response.get("result",{}).get("variantDbIds",[]):
                                    variantOutputLine = variants.get(variantDbId,["","",variantDbId,"","","","","",[]])
                                    variantOutputLines.append(variantOutputLine)
                            #initiate format
                            for entry in response.get("result",{}).get("dataMatrices",[]):
                                abbr = entry.get("dataMatrixAbbreviation",None)
                                if abbr and not abbr in variantOutputLine[8]:
                                    for variantOutputLine in variantOutputLines:
                                        variantOutputLine[8].append(abbr)
                                        for i in range(9,len(variantOutputLine)):
                                            variantOutputLine[i].append("")
                                    if not abbr in vcf_format_lines:
                                        vcf_format_lines[abbr] = ("##FORMAT=<ID=%s,Number=1,Type=%s,Description=\"%s\">" 
                                                % (abbr, str(entry.get("dataType","")), str(entry.get("dataMatrixName",""))))
                            #initiate output
                            nSamples = len(response.get("result",{}).get("callSetDbIds",[]))
                            for variantOutputLine in variantOutputLines:
                                variantOutputLine.extend([["."] * len(variantOutputLine[8])] * nSamples)
                            #populate output
                            for entry in response.get("result",{}).get("dataMatrices",[]):
                                abbr = entry.get("dataMatrixAbbreviation",None)
                                if abbr and abbr in variantOutputLine[8]:
                                    pos = variantOutputLine[8].index(abbr)
                                    for variantCounter in range(len(entry["dataMatrix"])):
                                        for sampleCounter in range(len(entry["dataMatrix"][variantCounter])):
                                            value = entry["dataMatrix"][variantCounter][sampleCounter]
                                            variantOutputLines[variantCounter][sampleCounter+sampleCounterOffset][pos] = value
                            sampleCounterOffset+=nSamples
                            callSetPage+=1
                        #output lines
                        for variantOutputLine in variantOutputLines:
                            for i in range(8,len(variantOutputLine)):
                                variantOutputLine[i] = ":".join(variantOutputLine[i])
                            buffer.write("%s\n" % "\t".join(variantOutputLine))
                        #next page
                        variantPage+=1
                    #finish output
                    for row in vcf_info_lines.values():
                        output.write("%s\n" % row)
                    for row in vcf_format_lines.values():
                        output.write("%s\n" % row)
                    buffer.seek(0)
                    shutil.copyfileobj(buffer,output)
            else:
                self._logger.error("No allelematrix in BrAPI endpoint")
        # except Exception as ex:
        #     self._logger.error(ex)
        finally:
            buffer.close()
        
