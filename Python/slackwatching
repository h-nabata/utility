# （Selenium 4.1.0 で動作確認済み）
# webdriver は使用しているChromeのバージョンに合わせて以下のページから取得する。
# https://chromedriver.chromium.org/downloads

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome import service as fs
import time
import sys
import os

# パラメータ
slackurl    = 'sample.slack.com'        # Slack URL
mailaddress = 'sample@gmail.com'        # 登録済みメールアドレス
mypassword  = 'abrakadabra'             # ログインパスワード
waittime    = 10.0                      # (sec.) 待機時間 (通信速度 & PCスペックに応じて調整)
logintime   = 10.0                      # (sec.) ログインを維持する時間 (適宜調整)
driverpath  = 'C:/Users/[user名]/Downloads/chromedriver_win32/chromedriver.exe' # webdriver が置かれている path を指定
workpath    = 'C:/Users/[user名]/Downloads/'                                    # 取得した画像などの保存先を指定
url         = 'https://app.slack.com/client/XXXYYYZZZ/AAABBBCCC'                # 目的地となる Slack WorkSpace のURL
recipient   = '山田'    # 相手(DMの受信者)のアカウント名
mymessage   = 'Yes My Darkness!'    # DMの内容
order       = 1  # 何番目に最新のメッセージを取得するか
dm_mode     = 0  # 0 -> 何もしない  //  1 -> メッセージを送信する
ss_mode     = 1  # 0 -> 何もしない  //  1 -> 最新のメッセージのスクリーンショットを保存  //  2 -> DM欄全体のスクリーンショットを保存


# フォルダの存在確認
if os.path.exists(workpath) == False:
    print('% Error: The "workpath" is invalid...')
    sys.exit()

# webdriverの存在確認
if os.path.exists(driverpath) == False:
    print('% Error: The "driverpath" is invalid...')
    sys.exit()

try:
    # 自動ログイン
    try:
        # Webページにアクセス
        chrome_service = fs.Service(executable_path=driverpath)
        driver  = webdriver.Chrome(service=chrome_service)
        driver.maximize_window()  #全画面表示にする (オプション)
        driver.get(url)
        driver.implicitly_wait(waittime)  # 遷移待ち

        # Slack WorkSpaceの選択
        ch_element = driver.find_element(By.ID, 'domain')
        ch_element.send_keys(slackurl)
        ch_element.send_keys(Keys.ENTER)
        driver.implicitly_wait(waittime)  # 遷移待ち

        # email, password の自動入力
        mail_element = driver.find_element(By.ID, 'email')
        mail_element.send_keys(mailaddress)
        driver.implicitly_wait(waittime)  # 待機 (オプション)
        password_element = driver.find_element(By.ID, 'password')
        password_element.send_keys(mypassword)
        password_element.send_keys(Keys.ENTER)
        driver.implicitly_wait(waittime)  # 遷移待ち (少し長め)

        # 自動ログインに成功
        print('% Successfully logged in:')
        url     = driver.current_url    # URLを取得する (オプション)
        title   = driver.title          # タイトルを取得する (オプション)
        print("% URL       :", url, )
        print("% title     :", title)

    except:
        # ログインに失敗した場合はセッションを終了する
        driver.quit()
        print("% Error: Your login has been failed...")
        sys.exit()

    # 送信先(recipient)のDMを開く
    recipient_list = driver.find_elements(By.CSS_SELECTOR, '[data-qa="channel_sidebar_name_' + recipient + '"]')
    if len(recipient_list) > 0:    # DMリストに名前が表示されている場合
        print("% DM is showing in the current list.")
        driver.find_element(By.CSS_SELECTOR, '[data-qa="channel_sidebar_name_' + recipient + '"]').click()
        driver.implicitly_wait(waittime)  # 遷移待ち
    else:    # DMリストに名前が表示されていない場合
        print("% DM is NOT showing in the current list.")
        searchbox_element = driver.find_element(By.CSS_SELECTOR, '[class="c-button-unstyled p-top_nav__search p-top_nav__search--windows-linux"]')
        searchbox_element.send_keys("@" + recipient)
        driver.implicitly_wait(waittime)  # 待機 (オプション)
        driver.find_element(By.CSS_SELECTOR, '[class="c-truncate c-truncate--break_words"]').click()
        driver.implicitly_wait(waittime)  # 待機 (オプション)

    # DMの要素をリストとして取得
    message_element_list = driver.find_elements(By.CSS_SELECTOR, '[class="c-message_kit__gutter"]')  # メッセージの要素を特定
    messagenum = -1*order
    target_massage = message_element_list[messagenum]  # 最新のメッセージをテキストで取得
    print("% sender    :", target_massage.find_element(By.CSS_SELECTOR, '[class="c-message__sender c-message_kit__sender"]').text) # 最新のメッセージの送信者
    print("% timestamp :", target_massage.find_element(By.CSS_SELECTOR, '[class="c-timestamp__label"]').text) # 最新のメッセージの送信時刻
    if target_massage.find_elements(By.CSS_SELECTOR, '[class="p-rich_text_section"]'):
        print("% message   :", target_massage.find_element(By.CSS_SELECTOR, '[class="p-rich_text_section"]').text, sep="\n") # 最新のメッセージの本文
    else:
        print("% (There is no text in the message.)")  # テキストを含まないメッセージの場合
        ss_mode = 1  # スクリーンショットの保存を強制する (オプション)

    # メッセージを送る (オプション)
    if dm_mode == 1:
        post_element = driver.find_element(By.CSS_SELECTOR, '[class="ql-editor ql-blank"]')
        post_element.click()
        driver.implicitly_wait(waittime)  # 待機 (オプション)
        post_element.send_keys(mymessage)  # メッセージ入力
        post_element.send_keys(Keys.ENTER)
        driver.implicitly_wait(waittime)  # 待機 (オプション)
        driver.find_element(By.CSS_SELECTOR, '[class="c-icon c-icon--paperplane-filled"]').click()
        print("% Your message was successfully sent.")

    # DMの画像を取得する (オプション)
    if ss_mode == 1:
        time.sleep(waittime)  # 強制的な待機 (オプション)
        png = message_element_list[messagenum].screenshot_as_png  # 最新のメッセージのみをpngで取得
        with open(workpath + 'image_from_slack.png', 'wb') as f:
            f.write(png)  # 画像の書き出し
            print("% The image of the latest message has been saved successfully.")
    elif ss_mode == 2:
        time.sleep(waittime)  # 強制的な待機 (オプション)
        dmfield_element = driver.find_elements(By.CSS_SELECTOR, '[class="p-workspace__primary_view"]')  # DM欄の要素を特定
        png = dmfield_element[messagenum].screenshot_as_png  #DM欄をpngで取得
        with open(workpath + 'image_from_slack.png', 'wb') as f:
            f.write(png)  # 画像の書き出し
            print("% The image of the messages being displayed has been saved successfully.")


    time.sleep(logintime)    # 任意時間ログイン状態を維持

finally:
    driver.quit()    # セッションの終了
    print("% Done.")
    sys.exit()
