charOfNuc = 'ACGT-'
massOfNuc = [135.128, 111.103, 151.128, 125.107, 100.000]
validStopCodon = ['TAA', 'TAG', 'TGA']
minimumCodon = 5
minimumPercentageOfCG = 30


def inputOutputFileNames():
  print('This program reports information about DNA')
  print('nucleotide sequences that may encode proteins.')
  inputFileName = input('Input file name? ')
  outputFileName = input('Output file name? ')

  return inputFileName, outputFileName


def fileRead(inputFileName):
  f = open(inputFileName, 'r')
  lines = f.readlines()
  f.close()
  return lines


def fileWrite(outputFileName, data):
  f = open(outputFileName, 'w')
  f.write(data)
  f.close()


def parseDnaValue(fileData):
  regionNamesArr = []
  nucleotidesArr = []

  for i in range(len(fileData)):
    if i % 2 == 0:
      regionNamesArr.append(fileData[i].strip())
    else:
      nucleotidesArr.append(fileData[i].strip().upper())

  return regionNamesArr, nucleotidesArr


def countNuc(nucleotides):
  nucCount = []
  for i in range(len(charOfNuc)):
    nucCount.append(0)

  for i in range(len(nucleotides)):
    if nucleotides[i] == charOfNuc[0]:
      nucCount[0] += 1
    elif nucleotides[i] == charOfNuc[1]:
      nucCount[1] += 1
    elif nucleotides[i] == charOfNuc[2]:
      nucCount[2] += 1
    elif nucleotides[i] == charOfNuc[3]:
      nucCount[3] += 1
    elif nucleotides[i] == charOfNuc[4]:
      nucCount[4] += 1

  return nucCount


def calMass(nucCount):
  totalMass = 0
  eachMass = []
  eachMassPercent = []

  for i in range(len(nucCount)):
    eachMass.append(massOfNuc[i] * nucCount[i])
    totalMass += massOfNuc[i] * nucCount[i]

  for mass in eachMass:
    eachMassPercent.append(round(mass / totalMass * 100, 1))

  return eachMassPercent, round(totalMass, 1)


def checkCodons(nucleotides):
  codonList = []

  nucleotides = nucleotides.replace('-', '')

  for i in range(0, len(nucleotides), 3):
    codonList.append(nucleotides[i:i + 3])

  return codonList


def chkIsProtein(codonList, eachMassPercent):
  stopCodon = codonList[len(codonList) - 1]

  if codonList[0] != 'ATG':
    return False
  if stopCodon not in validStopCodon:
    return False
  if len(codonList) < minimumCodon:
    return False
  if eachMassPercent[1] + eachMassPercent[2] < minimumPercentageOfCG:
    return False

  return True


def analysis(regionName, nucleotides):
  outputData = ''
  nucCount = countNuc(nucleotides)
  eachMassPercent, totalMass = calMass(nucCount)
  codonList = checkCodons(nucleotides)
  isProtein = chkIsProtein(codonList, eachMassPercent)

  outputData += 'Region Name: ' + str(regionName) + '\n'
  outputData += 'Nucleotides: ' + str(nucleotides) + '\n'
  outputData += 'Nuc. Counts: ' + str(nucCount[:-1]) + '\n'
  outputData += 'Total Mass%: ' + str(
    eachMassPercent[:-1]) + ' of ' + str(totalMass) + '\n'
  outputData += 'Codons List: [' + ', '.join(codonList) + ']' + '\n'
  outputData += 'Is Protein?: '
  outputData += 'YES' if isProtein else 'NO'
  outputData += '\n\n'

  return outputData


def main():
  inputFileName, outputFileName = inputOutputFileNames()
  regionNamesArr, nucleotidesArr = parseDnaValue(fileRead(inputFileName))
  outputData = ''

  for i in range(len(regionNamesArr)):
    outputData += analysis(regionNamesArr[i], nucleotidesArr[i])

  fileWrite(outputFileName, outputData)


main()